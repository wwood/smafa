use bird_tool_utils::clap_utils::add_clap_verbosity_flags;
use bird_tool_utils::clap_utils::set_log_level as set_log_level_bird_tool_utils;
use clap::*;

use std::env;

use smafa::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut app = build_cli();
    let matches = app.clone().get_matches();
    set_log_level(&matches, false);

    match matches.subcommand_name() {
        Some("query") => {
            let m = matches.subcommand_matches("query").unwrap();
            set_log_level(m, true);
            let db_root = m.get_one::<String>("database").unwrap();
            let query_fasta = m.get_one::<String>("query").unwrap();
            smafa::query(db_root, query_fasta)
        }
        Some("makedb") => {
            let m = matches.subcommand_matches("makedb").unwrap();
            set_log_level(m, true);
            let subject_fasta = m.get_one::<String>("input").unwrap();
            let database = m.get_one::<String>("database").unwrap();
            smafa::makedb(subject_fasta, database)
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
                .arg(arg!(-i --input <FILE> "Subject sequences to search against"))
                .arg(arg!(-d --database <FILE> "Output DB filename")),
        ))
        .subcommand(
            add_clap_verbosity_flags(
                Command::new("query")
                    .about("Search a database")
                    .arg(arg!(-d --database <FILE> "Output from makedb")),
            )
            .arg(arg!(-q --query <FILE> "Query sequences to search with")),
        )
}
