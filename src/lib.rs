#![cfg_attr(feature = "unstable", feature(test))]

pub mod fragment_clusterer;

extern crate bio;
#[macro_use]
extern crate log;

extern crate bincode;
extern crate serde;
#[macro_use]
extern crate serde_derive;
extern crate snap;

use std::str;

use bio::alphabets::dna;
use bio::data_structures::fmindex::*;
use bio::data_structures::suffix_array::{suffix_array, RawSuffixArray};
use bio::data_structures::bwt::*;
use bio::io::fasta;

use std::io::BufWriter;
use std::fs::File;
use std::time::Instant;
use std::collections::HashMap;
use std::collections::HashSet;
use std::collections::BTreeSet;

use std::process;

#[derive(Serialize, Deserialize)]
struct SaveableFMIndex {
    suffix_array: RawSuffixArray,
    bwt: BWT,
    less: Less,
    occ: Occ
}

#[derive(Serialize, Deserialize)]
struct SmafaDB {
    version: u8,
    fm_index: SaveableFMIndex
}

struct UnpackedDB {
    saveable_fm_index: SaveableFMIndex,
    text: Vec<u8>,
}

fn generate_unpacked_db(db_fasta: &str) -> UnpackedDB {
    let reader = fasta::Reader::new(File::open(db_fasta).unwrap());
    let mut sequence_length: usize = 0; // often 60
    let mut text = vec![];

    let mut now = Instant::now();
    let mut new_now;

    let mut i: u64 = 0;
    for record in reader.records() {
        let record = record.unwrap();
        if sequence_length == 0 {
            sequence_length = record.seq().len();
            debug!("Found first sequence length {}", sequence_length)
        } else if record.seq().len() != sequence_length {
            eprintln!("Found sequences of different lengths, all sequences must be of the same sequence.");
            process::exit(1);
        }
        text.extend(record.seq().iter().cloned());
        text.extend_from_slice(b"$");
        i+=1;
    }
    // Add a lexigraphically small character at the end to satisfy BWT
    // constraints. Do not use '\0' in case this is introduced into the bio-rust
    // BWT itself.
    text.extend_from_slice(b"#");
    new_now = Instant::now();
    info!("Read in {} sequences in {} seconds.", i, new_now.duration_since(now).as_secs());

    let alphabet = dna::n_alphabet();
    let before_index_generation_time = Instant::now();
    now = new_now;
    let sa = suffix_array(&text);
    new_now = Instant::now(); debug!("suffix array {:?}", new_now.duration_since(now)); now = new_now;
    let bwt1 = bwt(&text, &sa);
    new_now = Instant::now(); debug!("bwt {:?}", new_now.duration_since(now)); now = new_now;
    let less = less(&bwt1, &alphabet);
    new_now = Instant::now(); debug!("less {:?}", new_now.duration_since(now)); now = new_now;
    let occ = Occ::new(&bwt1, 3, &alphabet);
    new_now = Instant::now(); debug!("occ {:?}", new_now.duration_since(now));
    info!("Generated index in {} seconds.",
          Instant::now().duration_since(before_index_generation_time).as_secs());

    let saveable_fm_index = SaveableFMIndex {
        suffix_array: sa,
        bwt: bwt1,
        less: less,
        occ: occ
    };
    let unpacked = UnpackedDB {
        saveable_fm_index: saveable_fm_index,
        text: text
    };
    return unpacked;
}

pub fn makedb(db_root: &str, db_fasta: &str){
    let unpacked = generate_unpacked_db(db_fasta);
    let smafa_db = SmafaDB {
        version: 1,
        fm_index: unpacked.saveable_fm_index
    };
    let filename = db_root;
    let f = File::create(filename.clone()).unwrap();
    let writer = BufWriter::new(f);
    let mut snapper = snap::Writer::new(writer);
    debug!("Writing {}", filename);
    bincode::serialize_into(&mut snapper, &smafa_db, bincode::Infinite).unwrap();
    debug!("Finished writing.")
}


fn determine_sequence_length(text: &[u8]) -> usize {
    let mut count = 0;
    for (i, &item) in text.iter().enumerate() {
        if item == b'$' as u8 {
            count = i;
            break
        }
    }
    return count
}


pub fn query(db_root: &str, query_fasta: &str, max_divergence: u32,
             print_stream: &mut std::io::Write){
    let f = File::open(db_root).expect("file not found");
    let mut unsnapper = snap::Reader::new(f);
    debug!("Deserialising DB ..");
    let smafa_db: SmafaDB = bincode::deserialize_from(&mut unsnapper, bincode::Infinite).unwrap();
    debug!("Finished deserialising DB.");
    let fm_saveable = smafa_db.fm_index;

    let sa = fm_saveable.suffix_array;
    let bwt = fm_saveable.bwt;
    let less = fm_saveable.less;
    let occ = fm_saveable.occ;

    debug!("Inverting BWT ..");
    let text = bio::data_structures::bwt::invert_bwt(&bwt);
    debug!("Finished.");

    let mut num_hits: u64 = 0;
    {
        let query_printer = |hit: Hit| {
            num_hits += 1;
            print_stream.write_fmt(
                format_args!("{}\t{}\t{}\t{}\t{}\n",
                             hit.seq_id,
                             str::from_utf8(hit.query_sequence).unwrap(),
                             str::from_utf8(hit.hit_sequence).unwrap(),
                             hit.divergence,
                             hit.total)).unwrap();
        };
        query_with_everything(
            query_fasta, max_divergence, query_printer,
            &sa, &bwt, &less, &occ, &text);
    }
    info!("Printed {} hit(s).", num_hits);
}


struct Hit<'a, 'b, 'c> {
    seq_id: &'a str,
    query_sequence: &'b [u8],
    hit_sequence: &'c [u8],
    divergence: u32,
    total: u32
}


fn query_with_everything<T>(
    query_fasta: &str,
    max_divergence: u32,
    mut hit_processor: T,
    sa: &RawSuffixArray,
    bwt: &BWT,
    less: &Less,
    occ: &Occ,
    text: &Vec<u8>)
    where T: FnMut(Hit) {

    let fm = FMIndex::new(bwt, less, occ);

    let sequence_length = determine_sequence_length(text);
    debug!("Found sequence length {}", sequence_length);

    let reader = fasta::Reader::new(File::open(query_fasta).unwrap());
    for record in reader.records() {
        let seq = record.unwrap();
        let pattern = seq.seq();
        if pattern.len() != sequence_length {
            eprintln!("Sequence '{}' is not the same length as sequences in the database.", seq.id());
            process::exit(3);
        }
        let mut printed_seqs: HashSet<usize> = HashSet::new();
        for i in 0..((sequence_length-5) / 5) {
            let intervals = fm.backward_search(pattern[(5*i)..(5*i+5)].iter());
            let some = intervals.occ(sa);

            // The text is of length 60, plus the sentinal (if sequence_length
            // is 60 as original). So the sequence is located in the text at
            // text[pos/61 .. pos+60] I think.
            for pos in some {
                let sequence_index = pos / (sequence_length + 1);
                let start = sequence_index * (sequence_length + 1);
                if pos-start == 5*i && !printed_seqs.contains(&start) {
                    printed_seqs.insert(start);
                    let subject = &text[(start..(start+sequence_length))];
                    let mut divergence = 0;
                    let total = subject.len();
                    for (s,p) in subject.iter().zip(pattern) {
                        if s != p {
                            divergence = divergence + 1;
                        }
                    }
                    if divergence <= max_divergence {
                        let hit = Hit {
                            seq_id: seq.id(),
                            query_sequence: pattern,
                            hit_sequence: subject,
                            divergence: divergence,
                            total: total as u32
                        };
                        (hit_processor)(hit)
                    }
                }
            }
        }
    }
}

#[derive(Debug)]
struct BestClusterHit {
    divergence: u32,
    hit_sequence_index: u32
}
pub fn cluster(fasta_file: &str, max_divergence: u32, print_stream: &mut std::io::Write) {
    // For timing
    let mut now = Instant::now();
    let mut new_now;

    // query results are processed in order. We can recognize singletons as
    // sequences that only hit themselves.

    // make a database out of the sequences
    let db = generate_unpacked_db(fasta_file);
    new_now = Instant::now(); debug!("index generation time {:?}", new_now.duration_since(now)); now = new_now;

    let res = do_clustering(&db, &fasta_file, max_divergence);

    // print out the clusters
    new_now = Instant::now(); debug!("clustering search time {:?}", new_now.duration_since(now)); now = new_now;
    print_uc_file(
        &res.representative_sequence_ids,
        &res.sequence_index_to_best_hit,
        &res.sequence_names,
        determine_sequence_length(&db.text),
        print_stream);
    new_now = Instant::now(); debug!("printing time {:?}", new_now.duration_since(now));
}




struct ClusteringResults {
    representative_sequence_ids: BTreeSet<u32>,
    sequence_index_to_best_hit: HashMap<u32, BestClusterHit>,
    sequence_names: Vec<Box<String>>
}

fn do_clustering<'a>(db: &'a UnpackedDB, fasta_file: &'a str, max_divergence: u32) -> ClusteringResults {
    let sa = &db.saveable_fm_index.suffix_array;
    let bwt = &db.saveable_fm_index.bwt;
    let less = &db.saveable_fm_index.less;
    let occ = &db.saveable_fm_index.occ;
    let text = &db.text;

    let fm = FMIndex::new(bwt, less, occ);

    let sequence_length = determine_sequence_length(text);
    debug!("Found sequence length {}", sequence_length);

    let mut sequence_index_to_best_hit: HashMap<u32, BestClusterHit> = HashMap::new();
    let mut representative_sequence_indexes: BTreeSet<u32> = BTreeSet::new();
    let mut sequence_names = vec![];

    let reader = fasta::Reader::new(File::open(fasta_file).unwrap());
    let mut query_sequence_index: u32 = 0;
    for record in reader.records() {
        let seq = record.unwrap();
        let pattern = seq.seq();
        if pattern.len() != sequence_length {
            eprintln!("Sequence '{}' is not the same length as sequences in the database.", seq.id());
            process::exit(3);
        }
        sequence_names.push(Box::new(String::from(seq.id())));
        let mut printed_seqs: HashSet<usize> = HashSet::new();
        for i in 0..((sequence_length-5) / 5) {
            let intervals = fm.backward_search(pattern[(5*i)..(5*i+5)].iter());
            let some = intervals.occ(sa);

            // The text is of length 60, plus the sentinal (if sequence_length
            // is 60 as original). So the sequence is located in the text at
            // text[pos/61 .. pos+60] I think.
            for pos in some {
                let sequence_index = pos / (sequence_length + 1);
                if sequence_index <= query_sequence_index as usize {
                    let start = sequence_index * (sequence_length + 1);
                    if pos-start == 5*i &&
                        !printed_seqs.contains(&start) {
                            printed_seqs.insert(start);
                            let subject = &text[(start..(start+sequence_length))];
                            let mut divergence = 0;
                            for (s,p) in subject.iter().zip(pattern) {
                                if s != p {
                                    divergence = divergence + 1;
                                }
                            }
                            if divergence <= max_divergence &&
                                representative_sequence_indexes.contains(&(sequence_index as u32)) {
                                    // if it is better than the current hit, replace it
                                    if sequence_index_to_best_hit.contains_key(&query_sequence_index){
                                        // work out the true/false for the if statement here to
                                        // get around the borrow checker.
                                        let mut found_a_better_one = false;
                                        {
                                            let last_best = sequence_index_to_best_hit.get(&query_sequence_index).unwrap();
                                            if divergence < last_best.divergence ||
                                                (divergence == last_best.divergence &&
                                                 (sequence_index as u32) < last_best.hit_sequence_index) {
                                                    found_a_better_one = true;
                                                }
                                        }
                                        if found_a_better_one {
                                            sequence_index_to_best_hit.insert(
                                                query_sequence_index,
                                                BestClusterHit {
                                                    divergence: divergence,
                                                    hit_sequence_index: sequence_index as u32
                                                });
                                        };
                                    } else {
                                        // else if there is no hit in the hashmap, add it
                                        sequence_index_to_best_hit.insert(
                                            query_sequence_index,
                                            BestClusterHit {
                                                divergence: divergence,
                                                hit_sequence_index: sequence_index as u32
                                            });
                                    }
                                }
                        }
                }
            }
        }
        // Account for singletons
        if !sequence_index_to_best_hit.contains_key(&query_sequence_index) {
            representative_sequence_indexes.insert(query_sequence_index);
        }
        query_sequence_index += 1;
    }

    return ClusteringResults {
        representative_sequence_ids: representative_sequence_indexes,
        sequence_index_to_best_hit: sequence_index_to_best_hit,
        sequence_names: sequence_names
    }
}


fn print_uc_file(
    representative_sequence_ids: &BTreeSet<u32>,
    sequence_index_to_best_hit: &HashMap<u32, BestClusterHit>,
    sequence_names: &Vec<Box<String>>,
    alignment_length: usize,
    print_stream: &mut std::io::Write) {

    // Print out cluster centres
    let mut rep_id_to_cluster_id = HashMap::new();
    {
        // work out the number of entries in each cluster
        let mut num_in_each_cluster = HashMap::new();
        for (_, best_hit) in sequence_index_to_best_hit {
            let rep_id = best_hit.hit_sequence_index;
            if num_in_each_cluster.contains_key(&rep_id) {
                let new_num = num_in_each_cluster.get(&rep_id).unwrap() + 1;
                num_in_each_cluster.insert(rep_id, new_num);
            } else {
                num_in_each_cluster.insert(rep_id, 1);
            };
        };
        let mut cluster_id = 0;
        for rep_id in representative_sequence_ids {
            let num_in_cluster = match num_in_each_cluster.get(rep_id) {
                Some(num) => *num+1,
                None => 1
            };
            writeln!(print_stream, "C\t{}\t{}\t*\t*\t*\t*\t*\t{}\t*",
                     cluster_id,
                     num_in_cluster,
                     sequence_names[*rep_id as usize]).unwrap();
            rep_id_to_cluster_id.insert(*rep_id, cluster_id);
            cluster_id += 1;
        }
    }

    // Print out cluster memberships (H) ones
    for (member, best_hit) in sequence_index_to_best_hit {
        let rep_id = best_hit.hit_sequence_index;
        // Round to 1 decimal place.
        let perc_id = ((alignment_length as u32 -best_hit.divergence) as f32 /
                       alignment_length as f32 * 1000.0).round() / 10.0;
        writeln!(print_stream, "H\t{}\t{}\t{}\t*\t*\t*\t{}M\t{}\t{}",
                 rep_id_to_cluster_id.get(&rep_id).unwrap(),
                 alignment_length,
                 perc_id,
                 alignment_length,
                 sequence_names[*member as usize],
                 sequence_names[rep_id as usize]).unwrap();
    }
}

#[cfg(test)]
mod tests {
    extern crate tempfile;
    use super::*;
    use std::io::Read;

    #[test]
    fn makedb_and_query100seqs() {
        let tf_db: tempfile::NamedTempFile = tempfile::NamedTempFile::new().unwrap();
        makedb(tf_db.path().to_str().unwrap(), "test/data/4.08.ribosomal_protein_L3_rplC.100random.fna");
        let mut res = vec!();
        query(tf_db.path().to_str().unwrap(), "test/data/4.08.ribosomal_protein_L3_rplC.100random.fna", 5, &mut res);
        let mut expected: String = "".to_lowercase();
        File::open("test/data/4.08.ribosomal_protein_L3_rplC.100random.fnaVitself.expected").
            unwrap().read_to_string(&mut expected).unwrap();
        assert_eq!(expected, String::from_utf8(res).unwrap());
    }

    #[test]
    fn test_cluster_hello_world() {
        let fasta = "test/data/random2_plus_last_like_first.fna";
        let mut res = vec!();
        cluster(fasta, 1, &mut res);
        let mut expected: String = "".to_lowercase();
        File::open("test/data/random2_plus_last_like_first.fna.cluster-divergence1.uc").
            unwrap().read_to_string(&mut expected).unwrap();
        assert_eq!(expected, String::from_utf8(res).unwrap());
    }

    #[test]
    fn test_cluster_all_singletons() {
        let fasta = "test/data/random2_plus_last_like_first.fna";
        let mut res = vec!();
        cluster(fasta, 0, &mut res);
        let mut expected: String = "".to_lowercase();
        File::open("test/data/random2_plus_last_like_first.fna.cluster-divergence0.uc").
            unwrap().read_to_string(&mut expected).unwrap();
        assert_eq!(expected, String::from_utf8(res).unwrap());
    }

    #[test]
    fn test_cluster_many_seqs() {
        let fasta = "test/data/singlem_plot_test.fna";
        let mut res = vec!();
        cluster(fasta, 5, &mut res);
        let mut expected: String = "".to_lowercase();
        // expected string not manually verified except that the correct number
        // of outputs is observed.
        File::open("test/data/singlem_plot_test.fna.uc").
            unwrap().read_to_string(&mut expected).unwrap();
        let mut expected_sorted: Vec<&str> = expected.split("\n").collect();
        let observed = String::from_utf8(res).unwrap();
        let mut observed_sorted: Vec<&str> = observed.split("\n").collect();
        assert_eq!(expected_sorted.sort(), observed_sorted.sort());
    }
}

#[cfg(all(feature = "unstable", test))]
mod benches {
    extern crate test;
    use super::*;

    #[bench]
    fn bench_sra_euk_seqs_query(b: &mut test::Bencher) {
        let fasta = "test/data/sra.euks.4.11.ribosomal_protein_L10.head1000.fa";
        let db = generate_unpacked_db(fasta);
        b.iter(|| {
            let mut num_hits = 0;
            let query_printer = |_hit: Hit| {
                num_hits += 1;
            };
            query_with_everything(
                fasta, 5, query_printer,
                &db.saveable_fm_index.suffix_array,
                &db.saveable_fm_index.bwt,
                &db.saveable_fm_index.less,
                &db.saveable_fm_index.occ,
                &db.text);
        });
    }

    #[bench]
    fn bench_sra_euk_seqs_clustering(b: &mut test::Bencher) {
        let fasta = "test/data/sra.euks.4.11.ribosomal_protein_L10.head1000.fa";
        let db = generate_unpacked_db(fasta);
        b.iter(|| do_clustering(
            &db,
            fasta,
            5
        ));
    }

}
