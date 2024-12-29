use eyre::ContextCompat;
use iset::IntervalSet;
use itertools::Itertools;
use noodles::{
    bed::{
        self,
        record::{Builder, OptionalFields},
    },
    core::Position,
    fasta::{record::Definition, Writer},
};
use std::{fs::File, io::Write};

use crate::utils::{generate_random_seq_ranges, write_misassembly};

#[derive(Debug, PartialEq, Eq, Clone)]
pub struct BrokenSequence {
    // Start of broken sequence.
    start: usize,
    // End of broken sequence.
    end: usize,
}

impl From<BrokenSequence> for Builder<3> {
    fn from(value: BrokenSequence) -> Self {
        let start = value.start.clamp(1, usize::MAX);
        bed::Record::<3>::builder()
            .set_start_position(Position::new(start).unwrap())
            .set_end_position(Position::new(value.end).unwrap())
            .set_optional_fields(OptionalFields::from(vec!["Broken".to_string()]))
    }
}

pub fn generate_breaks<'a>(
    seq: &'a str,
    regions: &IntervalSet<Position>,
    number: usize,
    seed: Option<u64>,
) -> eyre::Result<(Vec<&'a str>, Vec<BrokenSequence>)> {
    // Number of seqs is equal to number of breaks + 1.
    // Start (-|-|-) Stop
    let mut seqs = Vec::with_capacity(number + 1);
    let mut breaks: Vec<BrokenSequence> = vec![];
    let seq_segments = generate_random_seq_ranges(seq.len(), regions, 1, number, seed, true)?
        .context("No sequence segments")?
        .collect_vec();
    let mut seq_iter = seq_segments.into_iter().peekable();

    // Add starting sequence before first break.
    if let Some((_, _, brange)) = seq_iter.peek() {
        let segment = &seq[..brange.start];
        breaks.push(BrokenSequence {
            start: 1,
            end: brange.start + 1,
        });
        seqs.push(segment);
    };

    while let Some((_, _, brange)) = seq_iter.next() {
        if let Some((_, _, next_brange)) = seq_iter.peek() {
            seqs.push(&seq[brange.start..next_brange.start]);
            breaks.push(BrokenSequence {
                start: brange.start + 1,
                end: next_brange.start + 1,
            })
        } else {
            seqs.push(&seq[brange.start..seq.len()]);
            breaks.push(BrokenSequence {
                start: brange.start + 1,
                end: seq.len() + 1,
            })
        }
    }

    Ok((seqs, breaks))
}

pub fn write_breaks<O, R, I>(
    record_name: &str,
    seq_region_pairs: (Vec<&str>, I),
    writer_fa: &mut Writer<O>,
    output_bed: &mut Option<bed::Writer<File>>,
) -> eyre::Result<()>
where
    O: Write,
    R: TryInto<Builder<3>> + Clone,
    I: IntoIterator<Item = R>,
{
    for (i, (seq, region)) in seq_region_pairs
        .0
        .into_iter()
        .zip(seq_region_pairs.1)
        .enumerate()
    {
        let new_definition = TryInto::<Builder<3>>::try_into(region.clone())
            .map(|b| b.set_reference_sequence_name(record_name).build())
            .map(|r| {
                if let Ok(r) = r {
                    Definition::new(
                        format!("{record_name}:{}-{}", r.start_position(), r.end_position()),
                        None,
                    )
                } else {
                    Definition::new(format!("{record_name}_ctg_{i}"), None)
                }
            })
            .unwrap_or(Definition::new(format!("{record_name}_ctg_{i}"), None));

        write_misassembly(
            seq.bytes().collect_vec(),
            std::iter::once(region),
            new_definition,
            writer_fa,
            output_bed.as_mut(),
        )?;
    }

    Ok(())
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_generate_breaks() {
        let seq = "AAAGGCCCGGCCCGGGGATTTTATTTTGGGCCGCCCAATTTAATTT";
        let regions = IntervalSet::from_iter(std::iter::once(
            Position::new(1).unwrap()..Position::new(seq.len()).unwrap(),
        ));

        let (seqs, breaks) = generate_breaks(seq, &regions, 3, Some(42)).unwrap();
        assert_eq!(
            seqs,
            ["AAAGGCCCGGCCCGGG", "GATTTTAT", "TTTGGGCCGCCCAATTTAAT", "TT"]
        );
        assert_eq!(
            breaks,
            [
                BrokenSequence { start: 1, end: 17 },
                BrokenSequence { start: 17, end: 25 },
                BrokenSequence { start: 25, end: 45 },
                BrokenSequence {
                    start: 45,
                    end: seq.len() + 1
                }
            ]
        );
        assert_eq!(seqs.join(""), seq)
    }
}
