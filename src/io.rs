use eyre::Context;
use itertools::Itertools;
use noodles::{
    bgzf::{self, IndexedReader},
    fasta::{
        self,
        record::{definition::Definition, Sequence},
        Record, Writer,
    },
};
use rust_lapper::{Interval, Lapper};
use std::{
    collections::HashMap,
    fs::File,
    io::{stdout, BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
};

use crate::{
    sequence::{SequenceSegment, SequenceType},
    utils::calculate_new_coords,
};

type Outfiles = (Box<dyn Write>, Option<BufWriter<File>>);

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
        .map(BufWriter::new);

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
    mut reader_bed: Option<BufReader<File>>,
) -> Option<HashMap<String, Lapper<usize, ()>>> {
    reader_bed.as_mut().map(|input_bed| {
        let mut regions: HashMap<String, Lapper<usize, ()>> = HashMap::new();

        for rec in input_bed.lines().map_while(Result::ok) {
            let Some((ctg_name, start, end)) =
                rec.trim().split('\t').collect_tuple::<(&str, &str, &str)>()
            else {
                log::error!("Invalid BED3 line: {rec}");
                continue;
            };
            let (Ok(start), Ok(stop)) = (start.parse::<usize>(), end.parse::<usize>()) else {
                log::error!("Invalid start or end coordinate on line {rec}");
                continue;
            };
            let region = Interval {
                start,
                stop,
                val: (),
            };
            regions
                .entry(ctg_name.to_string())
                .and_modify(|r| {
                    r.insert(region.clone());
                })
                .or_insert_with(|| Lapper::new(vec![region]));
        }
        regions
    })
}

pub fn write_new_fasta(
    ctg_name: &str,
    ctg_seq: &str,
    seqs: &[SequenceSegment],
    fa_writer: &mut Writer<Box<dyn Write>>,
) -> eyre::Result<()> {
    let mut num_breaks = 0;
    let mut records: HashMap<String, String> = HashMap::new();

    for seq in seqs {
        // Is a break, start new contig name.
        let ctg_name = if let SequenceType::Break = seq.typ {
            num_breaks += 1;
            continue;
        } else if num_breaks == 0 {
            ctg_name.to_owned()
        } else {
            format!("{ctg_name}_{num_breaks}")
        };
        let seq_slice = if let SequenceType::Good = seq.typ {
            &ctg_seq[seq.itv.start..seq.itv.stop]
        } else if let Some(misassembled_sequence) = &seq.itv.val {
            misassembled_sequence.as_str()
        } else {
            continue;
        };
        records
            .entry(ctg_name)
            .and_modify(|seq| seq.push_str(seq_slice))
            .or_insert_with(|| seq_slice.to_owned());
    }
    for (definition, sequence) in records.into_iter().sorted_by(|a, b| a.0.cmp(&b.0)) {
        fa_writer.write_record(&Record::new(
            Definition::new(definition, None),
            Sequence::from(sequence.into_bytes()),
        ))?;
    }
    Ok(())
}

pub fn write_misassembly_bed(
    ctg_name: &str,
    seqs: &[SequenceSegment],
    bed_writer: &mut BufWriter<File>,
) -> eyre::Result<()> {
    let new_coords = calculate_new_coords(seqs);
    for (seq, new_coords) in seqs.iter().zip(new_coords) {
        let (new_start, new_stop) = (new_coords.start, new_coords.stop);
        let color = seq.typ.as_color();
        writeln!(
            bed_writer,
            "{}\t{}\t{}\t{}\t{}\t.\t{}\t{}\t{}",
            ctg_name,
            seq.itv.start, // Original start
            seq.itv.stop,  // Original end
            seq.typ,       // Sequence type
            0,
            // No strand
            new_start, // Adjusted start
            new_stop,  // Adjusted end
            color,     // Red
        )?;
    }
    Ok(())
}

#[cfg(test)]
mod test {
    use std::{
        fs::File,
        io::{BufRead, BufReader, BufWriter, Write},
    };

    use itertools::Itertools;
    use noodles::fasta::{
        record::{Definition, Sequence},
        Reader, Record, Writer,
    };
    use rust_lapper::{Interval, Lapper};

    use crate::{
        cli::Misassembly,
        io::{write_misassembly_bed, write_new_fasta},
        misassembly::create_all_sequences,
        sequence::SequenceSegment,
    };

    fn example_seq() -> (&'static str, &'static str, Vec<SequenceSegment>) {
        (
            "seq_1",
            "ATTATTATTGCA",
            create_all_sequences(
                &Misassembly::Misjoin {
                    number: 2,
                    length: 4,
                },
                "ATTATTATTGCA",
                &Lapper::new(vec![Interval {
                    start: 1,
                    stop: 5,
                    val: (),
                }]),
                Some(12),
                true,
            )
            .unwrap(),
        )
    }

    #[test]
    fn test_write_fa() {
        // Create misjoin
        let (ctg_name, ctg_seq, seqs) = example_seq();

        // Write to tempfile
        let fname = "test/test.fa";
        {
            let fh = Box::new(BufWriter::new(File::create(fname).unwrap())) as Box<dyn Write>;
            let mut writer = Writer::new(fh);
            write_new_fasta(ctg_name, ctg_seq, &seqs, &mut writer).unwrap()
        };
        // Check that sequence is correctly deleted.
        {
            let mut fh = Reader::new(BufReader::new(File::open(fname).unwrap()));
            let record = fh.records().next().unwrap().unwrap();
            // Remove nts from 2-6 (TATT)
            assert_eq!(
                record,
                Record::new(
                    Definition::new(ctg_name, None),
                    Sequence::from("ATATTGCA".as_bytes().to_vec())
                )
            )
        };
        // Then remove.
        std::fs::remove_file(fname).unwrap();
    }

    #[test]
    fn test_write_bed() {
        let (ctg_name, _ctg_seq, seqs) = example_seq();

        let fname = "test/test.bed";
        {
            let mut writer = BufWriter::new(File::create(fname).unwrap());
            write_misassembly_bed(ctg_name, &seqs, &mut writer).unwrap();
        };
        {
            let fh = BufReader::new(File::open(fname).unwrap());
            let res = fh
                .lines()
                .map_while(Result::ok)
                .map(|l| l.trim().split('\t').map(|e| e.to_owned()).collect_vec())
                .collect_vec();
            assert_eq!(
                res,
                vec![
                    vec![
                        "seq_1".to_owned(),
                        "0".to_owned(),
                        "2".to_owned(),
                        "good".to_owned(),
                        "0".to_owned(),
                        ".".to_owned(),
                        "0".to_owned(),
                        "2".to_owned(),
                        "0,0,0".to_owned()
                    ],
                    vec![
                        "seq_1".to_owned(),
                        "2".to_owned(),
                        "3".to_owned(),
                        "misjoin".to_owned(),
                        "0".to_owned(),
                        ".".to_owned(),
                        "2".to_owned(),
                        "2".to_owned(),
                        "255,0,0".to_owned()
                    ],
                    vec![
                        "seq_1".to_owned(),
                        "3".to_owned(),
                        "6".to_owned(),
                        "misjoin".to_owned(),
                        "0".to_owned(),
                        ".".to_owned(),
                        "2".to_owned(),
                        "2".to_owned(),
                        "255,0,0".to_owned()
                    ],
                    vec![
                        "seq_1".to_owned(),
                        "6".to_owned(),
                        "12".to_owned(),
                        "good".to_owned(),
                        "0".to_owned(),
                        ".".to_owned(),
                        "2".to_owned(),
                        "8".to_owned(),
                        "0,0,0".to_owned()
                    ],
                ]
            )
        }
        std::fs::remove_file(fname).unwrap();
    }
}
