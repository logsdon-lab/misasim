use std::path::PathBuf;

use clap::{Parser, Subcommand};

#[derive(Parser)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

// TODO: This can be automatically generate with macros.
#[derive(Debug, PartialEq, Eq, Subcommand)]
pub enum Commands {
    /// Simulate a misjoin in a sequence.
    Misjoin {
        /// Number of misjoins to simulate.
        #[arg(short, long, default_value_t = 1)]
        number: usize,

        /// Max length of misjoin.
        #[arg(short, long, default_value_t = 5_000)]
        length: usize,

        /// Input sequence file.
        #[arg(short, long)]
        infile: PathBuf,

        /// Input bed file. Each region should map to a sequence from infile.
        #[arg(short = 'r', long)]
        inbedfile: Option<PathBuf>,

        /// Output sequence file.
        #[arg(short, long)]
        outfile: Option<PathBuf>,

        /// Output BED file with misassemblies.
        #[arg(short = 'b', long)]
        outbedfile: Option<PathBuf>,

        /// Seed to use for the random number generator.
        #[arg(short, long)]
        seed: Option<u64>,
    },

    /// Simulate a falsely duplicated sequence.
    FalseDuplication {
        /// Number of false duplications to simulate.
        #[arg(short, long, default_value_t = 1)]
        number: usize,

        /// Max length of sequence to duplicate.
        #[arg(short, long, default_value_t = 5_000)]
        length: usize,

        /// Maximum number of duplications for any single segment.
        #[arg(short, long, default_value_t = 3)]
        max_duplications: usize,

        /// Input sequence file.
        #[arg(short, long)]
        infile: PathBuf,

        /// Input bed file. Each region should map to a sequence from infile.
        #[arg(short = 'r')]
        inbedfile: Option<PathBuf>,

        /// Output sequence file.
        #[arg(short, long)]
        outfile: Option<PathBuf>,

        /// Output BED file with misassemblies.
        #[arg(short = 'b', long)]
        outbedfile: Option<PathBuf>,

        /// Seed to use for the random number generator.
        #[arg(short, long)]
        seed: Option<u64>,
    },

    /// Simulate a gap in a sequence.
    Gap {
        /// Number of gaps to simulate.
        #[arg(short, long, default_value_t = 1)]
        number: usize,

        /// Max length of gap simulate.
        #[arg(short, long, default_value_t = 5_000)]
        length: usize,

        /// Input sequence file.
        #[arg(short, long)]
        infile: PathBuf,

        /// Input bed file. Each region should map to a sequence from infile.
        #[arg(short = 'r')]
        inbedfile: Option<PathBuf>,

        /// Output sequence file.
        #[arg(short, long)]
        outfile: Option<PathBuf>,

        /// Output BED file with misassemblies.
        #[arg(short = 'b', long)]
        outbedfile: Option<PathBuf>,

        /// Seed to use for the random number generator.
        #[arg(short, long)]
        seed: Option<u64>,
    },

    /// Simulate a break in a sequence.
    Break {
        /// Number of breaks to simulate.
        #[arg(short, long, default_value_t = 1)]
        number: usize,

        /// Max length of break.
        #[arg(short, long, default_value_t = 5_000)]
        length: usize,

        /// Input sequence file.
        #[arg(short, long)]
        infile: PathBuf,

        /// Input bed file. Each region should map to a sequence from infile.
        #[arg(short = 'r')]
        inbedfile: Option<PathBuf>,

        /// Output sequence file.
        #[arg(short, long)]
        outfile: Option<PathBuf>,

        /// Output BED file with misassemblies.
        #[arg(short = 'b', long)]
        outbedfile: Option<PathBuf>,

        /// Seed to use for the random number generator.
        #[arg(short, long)]
        seed: Option<u64>,
    },
}
