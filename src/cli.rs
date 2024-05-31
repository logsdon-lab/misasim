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
    Collapse {
        /// Length in bp to collapse.
        /// Requires a repetitive sequence of this length or more in the input.
        #[arg(short, long)]
        length: usize
    },

    FalseDuplication {
        /// Length in bp to duplicate.
        #[arg(short, long)]
        length: usize
    }
}