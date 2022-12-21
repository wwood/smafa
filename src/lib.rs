use needletail::parse_fastx_file;
use serde::{Deserialize, Serialize};
use serde_json;
use std::io::{Read, Write};
use std::time::Instant;
use std::{error::Error, fs::File};

use log::{debug, info};

pub const AUTHOR_AND_EMAIL: &str =
    "Ben J. Woodcroft, Centre for Microbiome Research, School of Biomedical Sciences, Faculty of Health, Queensland University of Technology <benjwoodcroft near gmail.com>";

#[derive(Serialize, Deserialize, Debug)]
struct WindowSet {
    version: u32,
    windows: Vec<Vec<bool>>,
}

pub fn makedb(subject_fasta: &str, db_path: &str) -> Result<(), Box<dyn Error>> {
    // Iterate over lines, creating a vector of bools by one hot encoding
    // the input.

    // Open the query file as a fasta file.
    debug!("Opening subject fasta file: {}", subject_fasta);
    let mut subject_reader =
        parse_fastx_file(subject_fasta).expect("valid path/file of subject fasta");

    let mut encoded = Vec::new();
    info!("Encoding subject sequences ..");
    while let Some(record) = subject_reader.next() {
        let record = record.expect("valid record");
        let encoded1 = record
            .seq()
            .iter()
            .flat_map(|c| encode_single(*c))
            .collect::<Vec<_>>();
        encoded.push(encoded1);
    }

    let windows = WindowSet {
        version: 1,
        windows: encoded,
    };
    info!(
        "Encoding of {} sequences complete, writing db file {}",
        windows.windows.len(),
        db_path
    );

    // Encode
    let mut ferris_file = File::create(db_path)?;
    ferris_file.write_all(&postcard::to_allocvec(&windows).unwrap())?;
    info!("DB file written");
    Ok(())
}

// inline this function, performance affects untested, guessing it's better
#[inline(always)]
fn encode_single(c: u8) -> [bool; 5] {
    match c {
        b'A' => [true, false, false, false, false],
        b'C' => [false, true, false, false, false],
        b'G' => [false, false, true, false, false],
        b'T' => [false, false, false, true, false],
        b'N' | b'-' => [false, false, false, false, true],
        _ => {
            panic!("Invalid character in query sequence: {}", c)
        }
    }
}

fn get_hit_sequence(hit_bools: &Vec<bool>) -> String {
    let mut s = String::new();
    for j in 0..hit_bools.len() / 5 {
        let slice = &hit_bools[j * 5..(j + 1) * 5];
        s.push(match slice {
            [true, false, false, false, false] => 'A',
            [false, true, false, false, false] => 'C',
            [false, false, true, false, false] => 'G',
            [false, false, false, true, false] => 'T',
            [false, false, false, false, true] => 'N',
            _ => {
                panic!("Invalid character in query sequence: {:?}", slice)
            }
        });
    }
    s
}

pub fn query(
    db_path: &str,
    query_fasta: &str,
    max_divergence: Option<u32>,
    max_num_hits: Option<u32>,
) -> Result<(), Box<dyn Error>> {
    // Decode
    info!("Decoding db file {}", db_path);
    let start = Instant::now();
    let mut ferris_file = File::open(db_path)?;
    let mut buffer = Vec::new();
    ferris_file.read_to_end(&mut buffer)?;
    let windows: WindowSet = postcard::from_bytes(&buffer)?;

    // Open the query file as a fasta file.
    let mut query_reader = parse_fastx_file(query_fasta).expect("valid path/file of query fasta");

    // Iterate over the query file.
    info!("Querying ..");
    let mut query_number: u32 = 0;
    while let Some(record) = query_reader.next() {
        // encode a line from stdin as a vector of bools
        let query_vec = record
            .expect("Failed to parse query sequence")
            .seq()
            .iter()
            .flat_map(|c| encode_single(*c))
            .collect::<Vec<_>>();

        // Get the minimum distance between the query and each window using xor.
        let distances = windows
            .windows
            .iter()
            .map(|window| {
                window
                    .iter()
                    .zip(query_vec.iter())
                    .map(|(a, b)| a ^ b)
                    .filter(|x| *x)
                    .count()
            })
            .collect::<Vec<_>>();

        // Find the max_num_hits'th minimum distance.
        match max_num_hits {
            Some(max_num_hits) => {
                let mut min_distances = distances
                    .iter()
                    .enumerate()
                    .map(|(i, d)| (*d, i))
                    .collect::<Vec<_>>();
                // There might be a faster way of doing this using a priority
                // queue, but eh for now unless it really is slow.
                min_distances.sort();

                // If max num hits is greater than the number of windows, just print them all.
                let max_distance = match max_num_hits > min_distances.len() as u32 {
                    true => *distances.iter().max().unwrap(),
                    false => min_distances[(max_num_hits - 1) as usize].0,
                };

                // Print out the windows that qualify in order of increasing distance.
                for (distance, i) in min_distances.iter() {
                    if *distance <= max_distance
                        && (max_divergence.is_none()
                            || *distance / 2 <= max_divergence.unwrap() as usize)
                    {
                        let s = get_hit_sequence(&windows.windows[*i]);
                        println!("{}\t{}\t{}\t{}", query_number, i, distance / 2, s);
                    }
                }
            }
            None => {
                // Find the minimum distance.
                let min_distance = distances.iter().min().unwrap();
                debug!("Min distance: {}", min_distance);

                // Print the windows with the minimum distance.
                if max_divergence.is_none() || min_distance / 2 <= max_divergence.unwrap() as usize
                {
                    for (i, distance) in distances.iter().enumerate() {
                        if distance == min_distance {
                            let s = get_hit_sequence(&windows.windows[i]);
                            println!("{}\t{}\t{}\t{}", query_number, i, distance / 2, s);
                        }
                    }
                }
            }
        }

        query_number += 1;
    }

    info!(
        "Querying complete, took {} seconds",
        start.elapsed().as_secs()
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_makedb() {
        // Create a temporary directory to store the test DB file.
        let temp_dir = tempfile::tempdir().unwrap();
        let db_path = temp_dir.path().join("test.db");
        let db_path_str = db_path.to_str().unwrap();
        let subject_fasta = "tests/data/subjects.fa";

        // Call the makedb function with the test subject FASTA file and the path
        // to the test DB file.
        assert!(makedb(subject_fasta, &db_path_str).is_ok());

        // Check that the DB file exists.
        assert!(db_path.exists());

        // Open the DB file and decode it to a WindowSet struct.
        let mut ferris_file = File::open(&db_path).unwrap();
        let mut encoded = Vec::new();
        ferris_file.read_to_end(&mut encoded).unwrap();
        let windows = postcard::from_bytes::<WindowSet>(&encoded).unwrap();

        // Check that the WindowSet struct has the expected number of sequences.
        assert_eq!(windows.windows.len(), 5);

        // Check that the first sequence has the expected one-hot encoded values.
        let expected_encoded = vec![
            vec![true, false, false, false, false],
            vec![false, true, false, false, false],
            vec![false, false, true, false, false],
            vec![false, false, false, true, false],
            vec![false, false, false, false, true],
        ];
        assert_eq!(windows.windows, expected_encoded);
    }
}

// Derive IntoJson
#[derive(Serialize, Deserialize, Debug)]
struct CountResult {
    path: String,
    num_reads: usize,
    num_bases: usize,
}

pub fn count(paths: &Vec<&String>) -> Result<(), Box<dyn Error>> {
    let mut results = Vec::new();
    for path in paths {
        let mut reader = parse_fastx_file(path)?;
        let mut read_count = 0;
        let mut bases_count = 0;
        while let Some(record) = reader.next() {
            let record = record?;
            read_count += 1;
            bases_count += record.seq().len();
        }
        results.push(CountResult {
            path: path.to_string(),
            num_reads: read_count,
            num_bases: bases_count,
        });
    }
    // Print output in JSON format including input path
    println!("{}", serde_json::to_string(&results).unwrap());
    Ok(())
}
