use clap::{CommandFactory, Parser};
use collapse::generate_collapse;
use noodles::fasta::{reader::Records, Reader};
use std::{
    fs::File,
    io::{stdin, BufRead, BufReader, IsTerminal},
    path::PathBuf,
};

mod cli;
mod collapse;

use {
    cli::{Cli, Commands},
    collapse::find_all_repeats,
};

fn generate_misassemblies<B: BufRead>(
    records: Records<B>,
    command: cli::Commands,
) -> eyre::Result<()> {
    for record in records {
        let record = record?;
        let seq = std::str::from_utf8(record.sequence().as_ref())?;
        match command {
            cli::Commands::Misjoin { length, number } => {
                println!("Misjoin with length: {}", length);
            }
            cli::Commands::Collapse {
                length,
                num_repeats,
            } => {
                let repeats = find_all_repeats(seq, length);
                let seqs = generate_collapse(seq, &repeats, num_repeats);
            }
            cli::Commands::FalseDuplication { length, number } => {
                println!("False duplication with length: {}", length);
            }
            cli::Commands::Gap { length, number } => {
                println!("Gap with length: {}", length);
            }
        }
    }

    Ok(())
}

fn main() -> eyre::Result<()> {
    let (file, cmd) = if std::env::var("DEBUG").map_or(false, |v| v == "1" || v == "true") {
        let cmd = Commands::Collapse {
            length: 5,
            num_repeats: 1,
        };
        let file = PathBuf::from("test/data/test.fa");
        (file, cmd)
    } else {
        let cli = Cli::parse();
        let file = cli.infile;
        let cmd = cli.command;
        (file, cmd)
    };

    // https://rust-cli.github.io/book/in-depth/machine-communication.html
    if file == PathBuf::from("-") {
        if stdin().is_terminal() {
            Cli::command().print_help()?;
            std::process::exit(2);
        }

        let buf_reader = BufReader::new(stdin().lock());
        let mut fasta_reader = Reader::new(buf_reader);
        let results = generate_misassemblies(fasta_reader.records(), cmd);
    } else {
        let buf_reader = BufReader::new(File::open(&file).unwrap());
        let mut fasta_reader = Reader::new(buf_reader);
        let results = generate_misassemblies(fasta_reader.records(), cmd);
    }
    Ok(())
}
