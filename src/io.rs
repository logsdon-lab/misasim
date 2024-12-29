use eyre::Context;
use iset::IntervalSet;
use noodles::{
    bed,
    bgzf::{self, IndexedReader},
    core::Position,
    fasta,
};
use std::{
    collections::HashMap,
    fs::File,
    io::{stdout, BufReader, Write},
    path::{Path, PathBuf},
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

pub enum FastaReader {
    Bgzip(fasta::io::Reader<IndexedReader<File>>),
    Standard(fasta::io::Reader<BufReader<File>>),
}

pub struct Fasta {
    pub(crate) reader: FastaReader,
    pub(crate) index: fasta::fai::Index,
}

impl Fasta {
    pub fn new(infile: impl AsRef<Path>) -> eyre::Result<Self> {
        let (index, gzi) = Self::get_faidx(&infile)?;
        let fh = Self::read_fa(&infile, gzi.as_ref())?;
        Ok(Self { reader: fh, index })
    }

    pub fn lengths(&self) -> Vec<(String, u64)> {
        self.index
            .iter()
            .map(|rec| {
                (
                    String::from_utf8(rec.name().to_vec()).unwrap(),
                    rec.length(),
                )
            })
            .collect()
    }

    fn get_faidx(
        fa: &impl AsRef<Path>,
    ) -> eyre::Result<(fasta::fai::Index, Option<bgzf::gzi::Index>)> {
        // https://www.ginkgobioworks.com/2023/03/17/even-more-rapid-retrieval-from-very-large-files-with-rust/
        let fa_path = fa.as_ref().canonicalize()?;
        let is_bgzipped = fa_path.extension().and_then(|e| e.to_str()) == Some("gz");
        let fai_fname = fa_path.with_extension(if is_bgzipped { "gz.fai" } else { "fa.fai" });
        let fai = fasta::fai::read(fai_fname);
        if is_bgzipped {
            let index_reader = bgzf::indexed_reader::Builder::default()
                .build_from_path(fa)
                .with_context(|| format!("Failed to read gzi for {fa_path:?}"))?;
            let gzi = index_reader.index().clone();

            if let Ok(fai) = fai {
                log::debug!("Existing fai index found for {fa_path:?}");
                return Ok((fai, Some(gzi)));
            }
            log::debug!("No existing faidx for {fa_path:?}. Generating...");
            let mut records = Vec::new();
            let mut indexer = fasta::io::Indexer::new(index_reader);
            while let Some(record) = indexer.index_record()? {
                records.push(record);
            }

            Ok((fasta::fai::Index::from(records), Some(gzi)))
        } else {
            if let Ok(fai) = fai {
                return Ok((fai, None));
            }
            log::debug!("No existing faidx for {fa_path:?}. Generating...");
            Ok((fasta::index(fa)?, None))
        }
    }

    pub fn fetch(&mut self, ctg_name: &str, start: u32, stop: u32) -> eyre::Result<fasta::Record> {
        let start_pos = noodles::core::Position::new(start.clamp(1, u32::MAX) as usize).unwrap();
        let stop_pos = noodles::core::Position::new(stop.clamp(1, u32::MAX) as usize).unwrap();
        let region = noodles::core::Region::new(ctg_name, start_pos..=stop_pos);
        match &mut self.reader {
            FastaReader::Bgzip(reader) => Ok(reader.query(&self.index, &region)?),
            FastaReader::Standard(reader) => Ok(reader.query(&self.index, &region)?),
        }
    }

    fn read_fa(
        fa: &impl AsRef<Path>,
        fa_gzi: Option<&bgzf::gzi::Index>,
    ) -> eyre::Result<FastaReader> {
        let fa_file = std::fs::File::open(fa);
        if let Some(fa_gzi) = fa_gzi {
            Ok(FastaReader::Bgzip(
                fa_file
                    .map(|file| bgzf::IndexedReader::new(file, fa_gzi.to_vec()))
                    .map(fasta::io::Reader::new)?,
            ))
        } else {
            Ok(FastaReader::Standard(
                fa_file
                    .map(std::io::BufReader::new)
                    .map(fasta::io::Reader::new)?,
            ))
        }
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
