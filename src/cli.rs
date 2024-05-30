use std::path::PathBuf;

use clap::{Parser, Subcommand};


#[derive(Parser)]
pub struct Cli {
    /// Input sequence file.
    pub infile: PathBuf,

    /// Seed to use for the random number generator.
    #[arg(short, long)]
    pub seed: Option<u64>,

    #[command(subcommand)]
    pub command: Option<Commands>,
}


#[derive(Subcommand)]
pub enum Commands {
    /// Simulate a misjoin.
    Misjoin {
        /// Length in bp to remove.
        #[arg(short, long)]
        length: usize,
    },

    /// Simulate a collapse.
    /// Requires a repetitive sequence.
    Collapse,

    FalseDuplication {
        /// Length in bp to duplicate.
        #[arg(short, long)]
        length: usize
    }
}