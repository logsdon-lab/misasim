use crate::cli::Cli;
use clap::CommandFactory;
use iset::IntervalSet;
use log::info;
use noodles::{bed, bgzf, core::Position};
use std::{
    collections::HashMap,
    ffi::OsStr,
    fs::File,
    io::{stdin, stdout, BufRead, BufReader, IsTerminal, Write},
    path::PathBuf,
};

type Outfiles = (Box<dyn Write>, Option<bed::Writer<File>>);

pub fn get_outfile_writers(
    outfile: Option<PathBuf>,
    outbedfile: Option<PathBuf>,
) -> eyre::Result<Outfiles> {
    let output_fa: Box<dyn Write> = if let Some(outfile) = outfile {
        Box::new(File::create(outfile)?)
    } else {
        Box::new(stdout().lock())
    };
    let output_bed = outbedfile
        .and_then(|f| File::create(f).ok())
        .map(bed::Writer::new);

    Ok((output_fa, output_bed))
}

pub fn get_fa_reader(infile: PathBuf) -> eyre::Result<Box<dyn BufRead>> {
    if infile == PathBuf::from("-") {
        if stdin().is_terminal() {
            Cli::command().print_help()?;
            std::process::exit(2);
        }
        info!("Reading from stdin.");
        Ok(Box::new(BufReader::new(stdin().lock())) as Box<dyn BufRead>)
    } else if infile.extension() == Some(OsStr::new("gz")) {
        Ok(Box::new(bgzf::Reader::new(File::open(&infile)?)))
    } else {
        info!("Reading {infile:?}.");
        Ok(Box::new(BufReader::new(File::open(&infile)?)))
    }
}

pub fn get_regions(
    mut reader_bed: Option<bed::Reader<BufReader<File>>>,
) -> Option<HashMap<String, IntervalSet<Position>>> {
    reader_bed.as_mut().map(|input_bed| {
        let mut regions: HashMap<String, IntervalSet<Position>> = HashMap::new();
        for rec in input_bed.records::<3>().flatten() {
            let region = rec.start_position()..rec.end_position();
            regions
                .entry(rec.reference_sequence_name().to_string())
                .and_modify(|r| {
                    r.insert(region.clone());
                })
                .or_insert_with(|| {
                    let mut rs = IntervalSet::new();
                    rs.insert(region);
                    rs
                });
        }
        regions
    })
}
