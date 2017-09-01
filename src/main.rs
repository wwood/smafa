extern crate bio;
#[macro_use]
extern crate log;
use log::LogLevelFilter;
extern crate env_logger;
use env_logger::LogBuilder;

extern crate clap;
use clap::*;

use std::str;
use std::env;

use bio::alphabets::dna;
use bio::data_structures::fmindex::*;
use bio::data_structures::suffix_array::{suffix_array};
use bio::data_structures::bwt::*;
use bio::io::fasta;

use std::fs::File;
use std::time::Instant;
use std::collections::HashSet;
use std::process;

fn main(){
    let matches = App::new("smafa")
        .version("0.1.0")
        .author("Ben J. Woodcroft <benjwoodcroft near gmail.com>")
        .about("Read aligner for small pre-aligned sequences")
        .args_from_usage(
            "<DB_FASTA>          'Subject sequences to search against'
             <QUERY_FASTA>       'Query sequences to search with'
             -d, --divergence=[INTEGER] 'Maximum number of mismatches in reported hits [default: 5]'
             -v, --verbose       'Print extra debug logging information'
             -q, --quiet         'Unless there is an error, do not print logging information'")
    .get_matches();

    // Setup logging.
    let mut log_level = LogLevelFilter::Info;
    if matches.is_present("verbose") {
        log_level = LogLevelFilter::Debug;
    }
    if matches.is_present("quiet") {
        log_level = LogLevelFilter::Error;
    }
    let mut builder = LogBuilder::new();
    builder.filter(None, log_level);
    if env::var("RUST_LOG").is_ok() {
        builder.parse(&env::var("RUST_LOG").unwrap());
    }
    builder.init().unwrap();

    let db_fasta = matches.value_of("DB_FASTA").unwrap();
    let query_fasta = matches.value_of("QUERY_FASTA").unwrap();
    let max_divergence = value_t!(matches.value_of("divergence"), u32).unwrap_or(5);
    query(db_fasta, query_fasta, max_divergence);
}

fn query(db_fasta: &str, query_fasta: &str, max_divergence: u32){
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
                    if pos-start == 5*i {
                        printed_seqs.insert(start);
                        let subject = &text[(start..(start+60))];
                        let mut divergence = 0;
                        let mut total = 0;
                        for i in 0..60 {
                            if subject[i] != pattern[i] {
                                divergence = divergence + 1;
                                if divergence > max_divergence {
                                    break
                                }
                            }
                            total = total + 1;
                        }
                        if divergence <= max_divergence {
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
