use needletail::parse_fastx_file;
use serde::{Deserialize, Serialize};
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

pub fn query(db_path: &str, query_fasta: &str) -> Result<(), Box<dyn Error>> {
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

        // Find the minimum distance.
        let min_distance = distances.iter().min().unwrap();

        // Print the windows with the minimum distance.
        for (i, distance) in distances.iter().enumerate() {
            if distance == min_distance {
                let mut s = String::new();
                for j in 0..windows.windows[i].len() / 5 {
                    let slice = &windows.windows[i][j * 5..(j + 1) * 5];
                    s.push(match slice {
                        [true, false, false, false, false] => 'A',
                        [false, true, false, false, false] => 'C',
                        [false, false, true, false, false] => 'G',
                        [false, false, false, true, false] => 'T',
                        _ => 'N', // possibly a gap, but we don't know. eh.
                    });
                }
                println!("{}\t{}\t{}\t{}", query_number, i, distance, s);
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
