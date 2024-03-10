use needletail::parse_fastx_file;
use serde::{Deserialize, Serialize};

use std::io::{Read, Write};
use std::num::{NonZeroU8, NonZeroUsize};
use std::path::{Path, PathBuf};
use std::time::Instant;
use std::{error::Error, fs::File};

use log::{debug, info};

mod cluster;
pub use cluster::cluster;

pub const AUTHOR_AND_EMAIL: &str =
    "Ben J. Woodcroft, Centre for Microbiome Research, School of Biomedical Sciences, Faculty of Health, Queensland University of Technology <benjwoodcroft near gmail.com>";

pub const CURRENT_DB_VERSION: u32 = 2;

#[derive(Serialize, Deserialize, Debug, Clone)]
struct SeqEncoding(Vec<u64>);

#[derive(Clone)]
struct SeqEncodingLength {
    encoding: SeqEncoding,
    len: usize,
}

impl SeqEncodingLength {
    fn from_bytes(identifier: &[u8], seq: &[u8]) -> Self {
        // We encode chunks of 12 nucleotides to a u64
        let encoding = seq.chunks(12).enumerate().map(|(chunk_num, chunk)| {
            // The index i here is just to throw a good error message
            chunk.iter().enumerate().fold(0u64, |acc, (i, &byte)| {
                let b = encode_single(byte)
                    .unwrap_or_else(|| {
                        let seqname = String::from_utf8_lossy(identifier);
                        panic!(
                            "Byte {} cannot be interpreted as nucleotide, in sequence \"{}\" at position {}",
                            byte, seqname, 12 * chunk_num + i
                        )
                    })
                    .get();
                acc | ((b as u64) << (5 * i))
            })
        }).collect();
        Self {
            encoding: SeqEncoding(encoding),
            len: seq.len(),
        }
    }
}

#[derive(Serialize, Deserialize, Debug)]
struct WindowSet {
    version: u32,
    windows: Vec<SeqEncoding>,
    // None if there are no windows. Else, they must all be the same size
    len: Option<NonZeroUsize>,
}

impl WindowSet {
    fn new(version: u32) -> Self {
        WindowSet {
            version,
            windows: Vec::new(),
            len: None,
        }
    }

    fn get_distances(&self, seq: &SeqEncodingLength, distances: &mut [usize]) {
        if let Some(n) = self.len {
            if n.get() != seq.len {
                panic!(
                    "{}",
                    &format!("Cannot compute distances between seq of length {} and windows of lengths {}", seq.len, n.get())
                )
            }
        }
        for (window, distance) in self.windows.iter().zip(distances.iter_mut()) {
            *distance = window
                .0
                .iter()
                .zip(seq.encoding.0.iter())
                .map(|(a, b)| (a ^ b).count_ones() as usize)
                .sum::<usize>()
                / 2
        }
    }

    fn push_encoding(&mut self, encoding: SeqEncodingLength) {
        if let Some(n) = self.len {
            if n.get() != encoding.len {
                panic!(
                    "{}",
                    &format!(
                        "WindowSet seq length is {}, got a new sequence of length {}",
                        n, encoding.len
                    )
                )
            }
        } else {
            self.len = Some(
                encoding
                    .len
                    .try_into()
                    .expect("Cannot add empty sequence to WindowSet"),
            );
        }
        self.windows.push(encoding.encoding)
    }

    fn get_as_string(&self, index: usize) -> String {
        let uints = &self.windows[index].0;
        let v = (0..self.len.map(NonZeroUsize::get).unwrap_or(0))
            .map(|i| {
                let d = i / 12;
                let r = i % 12;
                let b = ((uints[d] >> (5 * r)) & 31) as u8;
                match b {
                    0b10000 => b'A',
                    0b01000 => b'C',
                    0b00100 => b'G',
                    0b00010 => b'T',
                    0b00001 => b'N',
                    _ => {
                        panic!("Invalid character in query sequence: {b}")
                    }
                }
            })
            .collect();
        // Safety: All the bytes above are ASCII, so it will never fail
        unsafe { String::from_utf8_unchecked(v) }
    }
}

pub fn makedb(subject_fasta: &Path, db_path: &Path) -> Result<(), Box<dyn Error>> {
    // Iterate over lines, creating a vector of u8, where the lowest 5 bits of the u8
    // are the input nucleotides, one-hot encoded

    // Open the query file as a fasta file.
    debug!("Opening subject fasta file: {:?}", subject_fasta);
    let mut subject_reader =
        parse_fastx_file(subject_fasta).expect("valid path/file of subject fasta");

    info!("Encoding subject sequences ..");
    let mut windows = WindowSet::new(CURRENT_DB_VERSION);
    while let Some(record) = subject_reader.next() {
        let record = record.expect("valid record");
        let encoded = SeqEncodingLength::from_bytes(record.id(), &record.seq());
        windows.push_encoding(encoded);
    }

    info!(
        "Encoding of {} sequences complete, writing db file {}",
        windows.windows.len(),
        db_path.to_string_lossy()
    );

    // Encode
    let mut ferris_file = File::create(db_path)?;
    ferris_file.write_all(&postcard::to_allocvec(&windows).unwrap())?;
    info!("DB file written");
    Ok(())
}

const fn create_lut() -> [u8; 256] {
    let mut lut = [0; 256];
    let mut i = 0;
    while i < 256 {
        let b = match i as u8 {
            b'A' | b'a' => 0b10000,
            b'C' | b'c' => 0b01000,
            b'G' | b'g' => 0b00100,
            b'T' | b't' | b'U' | b'u' => 0b00010,
            b'N' | b'W' | b'S' | b'M' | b'K' | b'R' | b'Y' | b'B' | b'D' | b'H' | b'V' | b'-'
            | b'n' | b'w' | b's' | b'm' | b'k' | b'r' | b'y' | b'b' | b'd' | b'h' | b'v' => 0b00001,
            _ => 0,
        };
        lut[i] = b;
        i += 1;
    }
    lut
}

const BYTE_LUT: [u8; 256] = create_lut();

// inline this function, performance affects untested, guessing it's better
#[inline(always)]
fn encode_single(c: u8) -> Option<NonZeroU8> {
    let lut: [u8; 256] = BYTE_LUT; // statically verify lut has 256 elements

    // Safety: We just verified it has indices 0-255, so a u8 can't be out of bounds
    let encoding = unsafe { *lut.get_unchecked(c as usize) };
    NonZeroU8::new(encoding)
}

pub fn query(
    db_path: &Path,
    query_fasta: &Path,
    max_divergence: Option<u32>,
    max_num_hits: Option<u32>,
    limit_per_sequence: Option<u32>,
) -> Result<(), Box<dyn Error>> {
    // Decode
    info!("Decoding db file {:?}", db_path);
    let start = Instant::now();
    let mut ferris_file = File::open(db_path)?;
    let mut buffer = Vec::new();
    ferris_file.read_to_end(&mut buffer)?;

    // Check that the version of the db file is the most recent. We do not
    // support backwards compatibility.
    let version: u32 = postcard::from_bytes(&buffer[0..4])?;
    if version != CURRENT_DB_VERSION {
        panic!("Unsupported db file version: {}. This version of smafa only works with version {} databases. The last version to support version 1 databases was v0.6.0.", version, CURRENT_DB_VERSION);
    }
    let windows: WindowSet = postcard::from_bytes(&buffer)?;

    // Open the query file as a fasta file.
    let mut query_reader = parse_fastx_file(query_fasta).expect("valid path/file of query fasta");

    // 1 is a special case, it is equivalent to None.
    let max_divergence_for_match = max_num_hits.filter(|&max_num_hits| max_num_hits != 1);

    // Pre-initialise the distances vector so don't have to continually reallocate.
    let mut distances = vec![0; windows.windows.len()];

    // Iterate over the query file.
    info!("Querying ..");
    let mut query_number: u32 = 0;
    while let Some(record) = query_reader.next() {
        // encode a line from stdin as a vector of bools
        let record = record.expect("Failed to parse query sequence");
        let query_vec = SeqEncodingLength::from_bytes(record.id(), &record.seq());

        // Get the minimum distance between the query and each window using xor.
        windows.get_distances(&query_vec, &mut distances);

        // Find the max_num_hits'th minimum distance.
        match max_divergence_for_match {
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
                let mut last_sequence: Option<(String, u32)> = None;
                let mut new_last_sequence: Option<(String, u32)>; // to get around borrow checker
                for (distance, i) in min_distances.iter() {
                    if *distance <= max_distance
                        && (max_divergence.is_none()
                            || *distance <= max_divergence.unwrap() as usize)
                    {
                        let s = windows.get_as_string(*i);
                        debug!("Found hit sequence {} at distance {}", s, distance);

                        if let Some(limit_per_sequence_unwrapped) = limit_per_sequence {
                            // limit per sequence
                            match &last_sequence {
                                Some((last_seq, last_seq_count)) => {
                                    if last_seq == &s {
                                        if last_seq_count >= &limit_per_sequence_unwrapped {
                                            continue;
                                        } else {
                                            new_last_sequence =
                                                Some((s.clone(), last_seq_count + 1));
                                        }
                                    } else {
                                        new_last_sequence = Some((s.clone(), 1));
                                    }
                                }
                                None => {
                                    new_last_sequence = Some((s.clone(), 1));
                                }
                            }
                            last_sequence = new_last_sequence;
                        }

                        // Print the window if we make it here.
                        println!("{}\t{}\t{}\t{}", query_number, i, distance, s);
                    }
                }
            }
            None => {
                // Find the minimum distance.
                let min_distance = distances.iter().min().unwrap();
                debug!("Min distance: {}", min_distance);

                if limit_per_sequence.is_some() {
                    panic!("limit_per_sequence is implemented unless max_num_hits > 1. It can be implemented by analogy, just haven't gotten around to it.");
                }

                // Print the windows with the minimum distance.
                if max_divergence.is_none() || *min_distance <= max_divergence.unwrap() as usize {
                    for (i, distance) in distances.iter().enumerate() {
                        if distance == min_distance {
                            let s = windows.get_as_string(i);
                            println!("{}\t{}\t{}\t{}", query_number, i, distance, s);
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
    use std::str::FromStr;

    use super::*;

    #[test]
    fn test_makedb() {
        // Create a temporary directory to store the test DB file.
        let temp_dir = tempfile::tempdir().unwrap();
        let db_path = temp_dir.path().join("test.db");
        let subject_fasta = std::path::PathBuf::from_str("tests/data/subjects.fa").unwrap();

        // Call the makedb function with the test subject FASTA file and the path
        // to the test DB file.
        assert!(makedb(&subject_fasta, &db_path).is_ok());

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
            vec![0b10000],
            vec![0b01000],
            vec![0b00100],
            vec![0b00010],
            vec![0b00001],
        ];
        for (i, j) in expected_encoded.iter().zip(windows.windows.iter()) {
            assert_eq!(i, &j.0)
        }
    }
}

// Derive IntoJson
#[derive(Serialize, Deserialize, Debug)]
struct CountResult {
    path: PathBuf,
    num_reads: usize,
    num_bases: usize,
}

pub fn count<T: Iterator<Item = P>, P: AsRef<Path>>(paths: T) -> Result<(), Box<dyn Error>> {
    let mut results = Vec::new();
    for path in paths {
        let mut reader = parse_fastx_file(&path)?;
        let mut read_count = 0;
        let mut bases_count = 0;
        while let Some(record) = reader.next() {
            let record = record?;
            read_count += 1;
            bases_count += record.seq().len();
        }
        results.push(CountResult {
            path: path.as_ref().to_owned(),
            num_reads: read_count,
            num_bases: bases_count,
        });
    }
    // Print output in JSON format including input path
    println!("{}", serde_json::to_string(&results).unwrap());
    Ok(())
}
