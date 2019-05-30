extern crate bio;

use std::str;
use std::fs::File;
use std::collections::HashMap;
use std::collections::HashSet;
use std;

use bio::io::fasta;
use bio::io::fasta::FastaRead;

pub fn cluster_by_fragment(input_fasta_path: &str, max_divergence: u8,
                           print_stream: &mut std::io::Write){
    let mut record = fasta::Record::new();
    let mut reader = fasta::Reader::new(
        File::open(input_fasta_path).expect("Failed to open query fasta"));
    let mut num_fragments: u8 = 0;
    let mut first = true;
    let mut sequence_length: u32 = 0;
    // Storage of each cluster rep fragment's actual sequence.
    // [cluster_id][fragment][position] => u8
    let mut fragment_strings: Vec<Vec<Vec<u8>>> = vec!();
    // fragment index to hash of str to list of cluster ids with that string
    let mut fragment_hash: Vec<HashMap<Vec<u8>, Vec<u32>>> = vec!();
    let mut next_cluster_id: u32 = 0;
    let mut seen_sequences: HashSet<Vec<u8>> = HashSet::new();
    let mut fragment_lengths: Vec<usize> = vec!();

    while reader.read(&mut record).is_ok() {
        let pattern: Vec<u8> = record.seq().to_vec();
        if pattern.len() == 0 {
            break; // finished all reads
        } else if !first && pattern.len() as u32 != sequence_length {
            panic!("Unexpectedly encountered sequence of length {}, expected {}",
                   pattern.len(), sequence_length);
        }
        if first {
            first = false;
            sequence_length = pattern.len() as u32;
            // determine the size of fragments - one more than the max_divergence
            num_fragments = max_divergence + 1;
            let fragment_size_floor = (sequence_length / num_fragments as u32) as u8;
            let mut num_extra_size = sequence_length - (fragment_size_floor * num_fragments) as u32;
            debug!("Divvying into {} fragments of size floor {}, and {} extras",
                   num_fragments, fragment_size_floor, num_extra_size);
            for _ in 0..num_fragments {
                if num_extra_size > 0 {
                    fragment_lengths.push((fragment_size_floor+1) as usize);
                    num_extra_size -= 1;
                } else {
                    fragment_lengths.push(fragment_size_floor as usize);
                }
            }
            debug!("Using fragment lengths {:?}", fragment_lengths);
            // setup fragment arrays
            for _ in 0..num_fragments {
                fragment_hash.push(HashMap::new())
            }
        }

        // Need to iterate over fragments from the current sequence
        let mut potential_hits: HashSet<u32> = HashSet::new();
        let mut current_fragments: Vec<Vec<u8>> = Vec::with_capacity(num_fragments as usize);
        {
            let mut fragment_counter: usize = 0;
            for l in &fragment_lengths {
                current_fragments.push(
                    pattern[fragment_counter .. (fragment_counter+l)]
                        .iter().map(|h| *h).collect());
                    fragment_counter += l;
            }
        }
        debug!("current fragments: {:?}", current_fragments);

        if seen_sequences.contains(&pattern) {
            debug!("Already seen this sequence, skipping");
            continue;
        }

        // for each sequence in the fasta file
        // for each fragment index
        for (i, subseq) in current_fragments.iter().enumerate() {
            // is sequence offset this fragment in the corresponding hash?
            // if so align add it to the list of sequences to compare with
            debug!("Testing to see if {:?} is in {:?}", subseq, fragment_hash[i]);
            match fragment_hash[i].get(subseq) {
                Some(cluster_ids) => {
                    for cluster_id in cluster_ids {
                        potential_hits.insert(*cluster_id);}},
                None => {}
            }
        }

        // if there is any sequences to compare to, compare to all of them.
        let mut best_cluster_id: Option<u32> = None;
        {
            let mut best_divergence: Option<u8> = None;
            for cluster_id in potential_hits.iter() {
                let mut divergence: u8 = 0;
                // if there are any within the max divergence, choose the one
                // with the lowest divergence, and in case of ties the lowest
                // cluster ID.
                for (i, subseq) in current_fragments.iter().enumerate() {
                    for (rep, query) in fragment_strings[*cluster_id as usize][i].iter()
                        .zip(subseq) {
                            if rep != query {
                                divergence += 1;
                            }
                    }
                }
                if divergence <= max_divergence &&
                    (best_divergence == None ||
                     divergence < best_divergence.unwrap() ||
                     (divergence == best_divergence.unwrap() && *cluster_id < best_cluster_id.unwrap())) {
                        debug!("Current best cluster for {:?} is now {}", current_fragments, *cluster_id);
                        best_divergence = Some(divergence);
                        best_cluster_id = Some(*cluster_id);
                }
            }
        }

        // if there are no hits within the max_divergence, add this sequence as
        // a cluster representative. Add each of its fragments to the hashes.
        match best_cluster_id {
            Some(cluster_id) => { // member of previous cluster
                print_stream.write_fmt(format_args!("{}\t", str::from_utf8(&pattern).unwrap())).unwrap();
                for i in 0..num_fragments {
                    write!(print_stream,
                        "{}", str::from_utf8(&fragment_strings[cluster_id as usize][i as usize]).unwrap())
                        .unwrap();
                }
                writeln!(print_stream).unwrap();
            },
            None => { // new cluster
                let cluster_id = next_cluster_id;
                next_cluster_id += 1;
                // Record this sequence
                let mut frags: Vec<Vec<u8>> = Vec::with_capacity(num_fragments as usize); // set capacity to save memory
                for (i, subseq) in current_fragments.into_iter().enumerate() {
                    frags.push(subseq.iter().map(|h| *h).collect());
                    // add to hash
                    if fragment_hash[i].contains_key(&subseq) {
                        fragment_hash[i].get_mut(&subseq).unwrap().push(cluster_id);
                    } else {
                        fragment_hash[i].insert(subseq, vec![cluster_id]);
                    }
                }
                fragment_strings.push(frags);
                let s = str::from_utf8(&pattern).unwrap();
                writeln!(print_stream, "{}\t{}", s, s).unwrap();
            }
        }

        seen_sequences.insert(pattern);
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_simple(){
        let mut stream = Cursor::new(Vec::new());
        cluster_by_fragment(
            "tests/data/cluster_dummy1.fna",
            1,
            &mut stream);
        assert_eq!(
            "ATGC\tATGC
ATGG\tATGC
AAAA\tAAAA
",
            str::from_utf8(stream.get_ref()).unwrap())
    }

    #[test]
    fn test_bug1(){
        let mut stream = Cursor::new(Vec::new());
        cluster_by_fragment(
            "tests/data/cluster_bug1.fna",
            2,
            &mut stream);
        assert_eq!(
            "ATGCAAAAA\tATGCAAAAA\n\
             ATAAAAAAA\tATGCAAAAA\n\
             TTAAAAAAA\tTTAAAAAAA\n",
            str::from_utf8(stream.get_ref()).unwrap())
    }

    #[test]
    fn test_best_hit_changes_bug(){
        // seq4 in the file shouldn't be reported otherwise there are two
        // sequences that are the same but are given different centroids.
        let mut stream = Cursor::new(Vec::new());
        cluster_by_fragment(
            "tests/data/cluster_best_hit_changes.fna",
            2,
            &mut stream);
        assert_eq!(
            "ATGCAAAAA\tATGCAAAAA\n\
             ATAAAAAAA\tATGCAAAAA\n\
             TTAAAAAAA\tTTAAAAAAA\n",
            str::from_utf8(stream.get_ref()).unwrap())
    }
}
