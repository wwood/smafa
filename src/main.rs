extern crate bio;
#[macro_use]
extern crate log;
use log::LogLevelFilter;
extern crate env_logger;
use env_logger::LogBuilder;

extern crate bincode;
extern crate serde;
#[macro_use]
extern crate serde_derive;
extern crate snap;

extern crate clap;
use clap::*;

use std::str;
use std::env;

use bio::alphabets::dna;
use bio::data_structures::fmindex::*;
use bio::data_structures::suffix_array::{suffix_array, RawSuffixArray};
use bio::data_structures::bwt::*;
use bio::io::fasta;

use std::io::BufWriter;
use std::fs::File;
use std::time::Instant;
use std::collections::HashSet;
use std::process;

fn main(){
    let mut app = build_cli();
    let matches = app.clone().get_matches();

    match matches.subcommand_name() {
        Some("query") => {
            let m = matches.subcommand_matches("query").unwrap();
            let db_root = m.value_of("DB").unwrap();
            let query_fasta = m.value_of("QUERY_FASTA").unwrap();
            let max_divergence = value_t!(m.value_of("divergence"), u32).unwrap_or(5);
            set_log_level(m);
            query(db_root, query_fasta, max_divergence);
        },
        Some("makedb") => {
            let m = matches.subcommand_matches("makedb").unwrap();
            let db_fasta = m.value_of("DB_FASTA").unwrap();
            let db_root = m.value_of("DB").unwrap();
            set_log_level(m);
            makedb(db_root, db_fasta);
        },
        _ => {
            app.print_help().unwrap();
            println!();
        }
    }
}

fn set_log_level(matches: &clap::ArgMatches) {
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
}

fn build_cli() -> App<'static, 'static> {
    let makedb_args: &'static str = "<DB_FASTA>  'Subject sequences to search against'
                       <DB>        'Output DB filename root'

                      -v, --verbose       'Print extra debug logging information'
                      -q, --quiet         'Unless there is an error, do not print logging information'";
    let query_args: &'static str = "<DB>           'Output from makedb'
                      <QUERY_FASTA> 'Query sequences to search with'
                      -d, --divergence=[INTEGER] 'Maximum number of mismatches in reported hits [default: 5]'

                      -v, --verbose       'Print extra debug logging information'
                      -q, --quiet         'Unless there is an error, do not print logging information'";

    return App::new("smafa")
        .version("0.1.0")
        .author("Ben J. Woodcroft <benjwoodcroft near gmail.com>")
        .about("Read aligner for small pre-aligned sequences")
        .args_from_usage("-v, --verbose       'Print extra debug logging information'
             -q, --quiet         'Unless there is an error, do not print logging information'")
        .subcommand(
            SubCommand::with_name("makedb")
                .about("Generate a searchable database")
                .args_from_usage(&makedb_args))
        .subcommand(
            SubCommand::with_name("query")
                .about("Search a database")
                .args_from_usage(&query_args));
}

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

fn makedb(db_root: &str, db_fasta: &str){
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
    text.extend_from_slice(b"\0");
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
    let smafa_db = SmafaDB {
        version: 1,
        fm_index: saveable_fm_index,
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


fn query(db_root: &str, query_fasta: &str, max_divergence: u32){
    let f = File::open(db_root).expect("file not found");
    let mut unsnapper = snap::Reader::new(f);
    let smafa_db: SmafaDB = bincode::deserialize_from(&mut unsnapper, bincode::Infinite).unwrap();
    let fm_saveable = smafa_db.fm_index;

    let sa = fm_saveable.suffix_array;
    let bwt = fm_saveable.bwt;
    let less = fm_saveable.less;
    let occ = fm_saveable.occ;

    let fm = FMIndex::new(&bwt, &less, &occ);
    debug!("Inverting BWT ..");
    let text = bio::data_structures::bwt::invert_bwt(&bwt);
    debug!("Finished.");

    let sequence_length = determine_sequence_length(&text);
    debug!("Found sequence length {}", sequence_length);


    let reader = fasta::Reader::new(File::open(query_fasta).unwrap());
    let mut num_hits: u64 = 0;
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
            let some = intervals.occ(&sa);

            // The text is of length 60, plus the sentinal (if sequence_length
            // is 60 as original). So the sequence is located in the text at
            // text[pos/61 .. pos+60] I think.
            for pos in some {
                let sequence_index = pos / (sequence_length + 1);
                let start = sequence_index * (sequence_length + 1);
                if !printed_seqs.contains(&start) {
                    if pos-start == 5*i {
                        printed_seqs.insert(start);
                        let subject = &text[(start..(start+sequence_length))];
                        let mut divergence = 0;
                        let mut total = 0;
                        for i in 0..sequence_length {
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
