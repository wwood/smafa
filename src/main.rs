extern crate smafa;

extern crate log;
extern crate env_logger;
use log::LevelFilter;
use env_logger::Builder;

extern crate clap;
use clap::*;

use std::env;

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
            let k = value_t!(m.value_of("kmer-length"), usize).unwrap_or(5);
            smafa::query(db_root, query_fasta, max_divergence, &mut std::io::stdout(), k);
        },
        Some("makedb") => {
            let m = matches.subcommand_matches("makedb").unwrap();
            let db_fasta = m.value_of("DB_FASTA").unwrap();
            let db_root = m.value_of("DB").unwrap();
            set_log_level(m);
            smafa::makedb(
                db_root,
                db_fasta,
                match m.is_present("amino-acid") {
                    true => smafa::DatabaseType::Translated,
                    false => smafa::DatabaseType::DNA
                });
        },
        Some("cluster") => {
            let m = matches.subcommand_matches("cluster").unwrap();
            let query_fasta = m.value_of("FASTA").unwrap();
            let max_divergence = value_t!(m.value_of("divergence"), u32).unwrap_or(5);
            set_log_level(m);
            if m.is_present("fragment-method") {
                smafa::fragment_clusterer::cluster_by_fragment(
                    query_fasta, max_divergence as u8, &mut std::io::stdout());
            } else {
                smafa::cluster(
                    query_fasta,
                    max_divergence,
                    &mut std::io::stdout(),
                    match m.is_present("amino-acid") {
                        true => smafa::DatabaseType::Translated,
                        false => smafa::DatabaseType::DNA
                    });
            }
        },
        _ => {
            app.print_help().unwrap();
            println!();
        }
    }
}

fn set_log_level(matches: &clap::ArgMatches) {
    let mut log_level = LevelFilter::Info;
    if matches.is_present("verbose") {
        log_level = LevelFilter::Debug;
    }
    if matches.is_present("quiet") {
        log_level = LevelFilter::Error;
    }
    let mut builder = Builder::new();
    builder.filter_level(log_level);
    if env::var("RUST_LOG").is_ok() {
        builder.parse_filters(&env::var("RUST_LOG").unwrap());
    }
    builder.try_init().unwrap();
}

fn build_cli() -> App<'static, 'static> {
    let makedb_args: &'static str = "<DB_FASTA>  'Subject sequences to search against'
                       <DB>        'Output DB filename root'
                       --amino-acid  'Sequences are amino acid [default: nucleotide]'

                      -v, --verbose       'Print extra debug logging information'
                      -q, --quiet         'Unless there is an error, do not print logging information'";
    let query_args: &'static str = "<DB>           'Output from makedb'
                      <QUERY_FASTA> 'Query sequences to search with'
                      -d, --divergence=[INTEGER] 'Maximum number of mismatches in reported hits [default: 5]'
                      -k, --kmer-length=[INTEGER]   'Length of kmer to query with [default 5]'

                      -v, --verbose       'Print extra debug logging information'
                      -q, --quiet         'Unless there is an error, do not print logging information'";
    let cluster_args: &'static str = "<FASTA> 'Sequences to cluster'
                      -d, --divergence=[INTEGER] 'Maximum number of mismatches in reported hits [default: 5]'
                       --amino-acid  'Sequences are amino acid [default: nucleotide]'

                      -v, --verbose       'Print extra debug logging information'
                      -q, --quiet         'Unless there is an error, do not print logging information'";

    return App::new("smafa")
        .version(env!("CARGO_PKG_VERSION"))
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
                .args_from_usage(&query_args))
        .subcommand(
            SubCommand::with_name("cluster")
                .about("Cluster sequences greedily, preferring sequences towards front of file")
                .arg(Arg::with_name("fragment-method")
                     .long("fragment-method")
                     .help("Use the 'fragment' method for clustering"))
                .args_from_usage(&cluster_args));
}

