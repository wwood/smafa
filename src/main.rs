extern crate bio;
#[macro_use]
extern crate log;
extern crate env_logger;
extern crate bincode;
extern crate serde;
extern crate serde_json;

use std::str;

use bio::alphabets::dna;
use bio::data_structures::fmindex::*;
use bio::data_structures::suffix_array::{suffix_array};
use bio::data_structures::bwt::*;
use bio::io::fasta;

use std::env;
use std::fs::File;
use std::time::Instant;
use std::collections::HashSet;
use std::process;

fn main(){
    let args: Vec<String> = env::args().collect();
    let _ = env_logger::init().unwrap();

    match args.len() {
        4 => {
            match &args[1] as &str {
                "query" => query(&args[2], &args[3]),
                _ => {
                    print_help();
                    eprintln!("Unexpected argument '{}'.", &args[1]);
                    process::exit(1);
                }
            }
        },
        _ => {
            print_help();
            eprintln!("Unexpected number of arguments.");
            process::exit(1);
        }
    }
}

fn print_help(){
    println!("Usage:");
    println!("");
    println!(" smafa query <reference_fasta> <query_fasta>");
    println!();
}

fn query(db_fasta: &str, query_fasta: &str){
    let mut text = vec![];

    let reader = fasta::Reader::new(File::open(db_fasta).unwrap());

    let mut now = Instant::now();
    let mut new_now;

    {
        let mut i: u64 = 0;
        for record in reader.records() {
            let record = record.unwrap();
            if record.seq().len() != 60 {
                eprintln!("Currently only sequences of exactly length 60 are permitted.");
                process::exit(1);
            }
            text.extend(record.seq().iter().cloned());
            text.extend_from_slice(b"$");
            i+=1;
        }
        info!("Read in {} sequences.", i)
    }

    let alphabet = dna::n_alphabet();
    let before_index_generation_time = Instant::now();
    new_now = Instant::now(); debug!("reading {:?}", new_now.duration_since(now)); now = new_now;
    let sa = suffix_array(&text);
    new_now = Instant::now(); debug!("suffix array {:?}", new_now.duration_since(now)); now = new_now;
    let bwt = bwt(&text, &sa);
    new_now = Instant::now(); debug!("bwt {:?}", new_now.duration_since(now)); now = new_now;
    let less = less(&bwt, &alphabet);
    new_now = Instant::now(); debug!("less {:?}", new_now.duration_since(now)); now = new_now;
    let occ = Occ::new(&bwt, 3, &alphabet);
    new_now = Instant::now(); debug!("occ {:?}", new_now.duration_since(now)); now = new_now;
    let fm = FMIndex::new(&bwt, &less, &occ);
    new_now = Instant::now(); debug!("fmindex {:?}", new_now.duration_since(now));
    info!("Generated index in {:?}", Instant::now().duration_since(before_index_generation_time));

    let reader = fasta::Reader::new(File::open(query_fasta).unwrap());
    let mut num_hits: u64 = 0;
    for record in reader.records() {
        let seq = record.unwrap();
        let pattern = seq.seq();
        let mut printed_seqs: HashSet<usize> = HashSet::new();
        for i in 0..11 {
            let intervals = fm.backward_search(pattern[(5*i)..(5*i+5)].iter());
            let some = intervals.occ(&sa);

            // The text if of length 60, plus the sentinal.
            // So the sequence is located in the text at text[pos/61 .. pos+60] I think.
            for pos in some {
                let start2 = pos / 61;
                let start = start2 * 61;
                if !printed_seqs.contains(&start) {
                    printed_seqs.insert(start);
                    if pos-start == 0 {
                        let subject = &text[(start..(start+60))];
                        let mut divergence = 0;
                        let mut total = 0;
                        for i in 0..60 {
                            if subject[i] != pattern[i] {
                                divergence = divergence + 1;
                                if divergence > 5 {
                                    break
                                }
                            }
                            total = total + 1;
                        }
                        if divergence <= 5 {
                            num_hits = num_hits + 1;
                            println!("{}\t{}\t{}\t{}\t{}",
                                     seq.id(),
                                     str::from_utf8(pattern).unwrap(),
                                     str::from_utf8(subject).unwrap(),
                                     divergence, total)
                        }
                    }
                }
            }
        }
    }
    info!("Printed {} hit(s).", num_hits)
}
