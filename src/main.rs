use bird_tool_utils::clap_utils::add_clap_verbosity_flags;
use bird_tool_utils::clap_utils::set_log_level as set_log_level_bird_tool_utils;
use clap::*;

use std::env;
use std::path::PathBuf;

use smafa::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut app = build_cli();
    let matches = app.clone().get_matches();
    set_log_level(&matches, false);

    match matches.subcommand_name() {
        Some("query") => {
            let m = matches.subcommand_matches("query").unwrap();
            set_log_level(m, true);
            let db_root = m.get_one::<PathBuf>("database").unwrap();
            let query_fasta = m.get_one::<PathBuf>("query").unwrap();
            let max_divergence = m.get_one::<u32>("max-divergence");
            let max_num_hits = m.get_one::<u32>("max-num-hits");
            let limit_per_sequence = m.get_one::<u32>("limit-per-sequence");
            smafa::query(
                db_root,
                query_fasta,
                max_divergence.copied(),
                max_num_hits.copied(),
                limit_per_sequence.copied(),
            )
        }
        Some("makedb") => {
            let m = matches.subcommand_matches("makedb").unwrap();
            set_log_level(m, true);
            let subject_fasta = m.get_one::<PathBuf>("input").unwrap();
            let database = m.get_one::<PathBuf>("database").unwrap();
            smafa::makedb(subject_fasta, database)
        }
        Some("cluster") => {
            let m = matches.subcommand_matches("cluster").unwrap();
            set_log_level(m, true);
            let input_fasta = m.get_one::<PathBuf>("input").unwrap();
            let max_divergence = m.get_one::<u32>("max-divergence").unwrap();
            smafa::cluster(input_fasta, *max_divergence, &mut std::io::stdout())
        }
        Some("count") => {
            let m = matches.subcommand_matches("count").unwrap();
            set_log_level(m, true);
            let fastx_files = m.get_many::<PathBuf>("input").unwrap().collect::<Vec<_>>();
            smafa::count(fastx_files.iter())
        }
        _ => {
            app.print_help().unwrap();
            println!();
            Ok(())
        }
    }
}

fn set_log_level(matches: &clap::ArgMatches, is_last: bool) {
    set_log_level_bird_tool_utils(matches, is_last, "Smafa", crate_version!());
}

fn build_cli() -> Command {
    command!()
        .author(crate::AUTHOR_AND_EMAIL)
        .arg(arg!(-v --verbose "Print extra debug logging information"))
        .arg(arg!(-q --quiet "Unless there is an error, do not print logging information"))
        .subcommand(add_clap_verbosity_flags(
            Command::new("makedb")
                .about("Generate a searchable database")
                .arg(arg!(-i --input <FILE> "Subject sequences to search against").value_parser(value_parser!(PathBuf)))
                .arg(arg!(-d --database <FILE> "Output DB filename").value_parser(value_parser!(PathBuf)))
        ))
        .subcommand(add_clap_verbosity_flags(
            Command::new("query")
                .about("Search a database. See query --help for more information about output format.")
                .long_about("This command searches a database for query sequences. The database must be generated with the `makedb` command. The query sequences can be in FASTA or FASTQ format. The output is a tab-separated file with the following columns:\n\
                \n\
                1. Query sequence number (0-indexed)\n\
                2. Subject sequence number (0-indexed)\n\
                3. Divergence (number of nucleotides different between the two sequences\n\
                4. Subject sequence (with dashes and degenerate base symbols shown as Ns)")
                .arg(arg!(-d --database <FILE> "Output from makedb [required]").value_parser(value_parser!(PathBuf)))
                .arg(arg!(-q --query <FILE> "Query sequences to search with in FASTX format [required]").value_parser(value_parser!(PathBuf)))
                .arg(
                    arg!( --"max-divergence" <INT> "Maximum divergence to report hits for, for each sequence [default: not used]")
                        .value_parser(value_parser!(u32)),
                )
                .arg(
                    arg!( --"max-num-hits" <INT> "Maximum number of hits to report [default: 1]")
                        .value_parser(value_parser!(u32)),
                )
                .arg(
                    arg!( --"limit-per-sequence" <INT> "Maximum number of hits to report per sequence. Requires --max-num-hits > 1 for now. [default: not used]")
                        .value_parser(value_parser!(u32)),
                ),
        ))
        .subcommand(add_clap_verbosity_flags(
            Command::new("cluster")
                .about("Cluster sequences by similarity")
                .arg(arg!(-i --input <FILE> "FASTA file to cluster").value_parser(value_parser!(PathBuf)))
                .arg(
                    arg!(-d --"max-divergence" <INT> "Maximum divergence to report hits for, for each sequence [default: not used]")
                        .value_parser(value_parser!(u32)),
                ),
        ))
        .subcommand(add_clap_verbosity_flags(
            Command::new("count")
                .about("Print the number of reads/bases in a possibly gzipped FASTX file")
                .arg(arg!(-i --input <FILE> "FASTQ file to count").value_parser(value_parser!(PathBuf))
                // allow multipl
                .num_args(0..)
            ),
        ))
}
