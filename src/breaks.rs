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

use crate::utils::{get_sequence_segments, write_misassembly};

#[derive(Debug, PartialEq, Eq)]
pub struct SequenceBreak(pub usize);

impl From<SequenceBreak> for Builder<3> {
    fn from(value: SequenceBreak) -> Self {
        let pos = value.0.clamp(1, usize::MAX);
        bed::Record::<3>::builder()
            .set_start_position(Position::new(pos).unwrap())
            .set_end_position(Position::new(pos).unwrap())
            .set_optional_fields(OptionalFields::from(vec!["Break".to_string()]))
    }
}

pub fn generate_breaks(
    seq: &str,
    number: usize,
    seed: Option<u64>,
) -> eyre::Result<(Vec<&str>, Vec<Option<SequenceBreak>>)> {
    // Number of seqs is equal to number of breaks + 1.
    // Start (-|-|-) Stop
    let mut seqs = Vec::with_capacity(number + 1);
    let mut breaks: Vec<Option<SequenceBreak>> = vec![];
    let seq_segments = get_sequence_segments(seq, number, seed)
        .unwrap()
        .collect_vec();

    // Add starting sequence before first break.
    if let Some((start, _, _)) = seq_segments.first() {
        let segment = &seq[..*start];
        seqs.push(segment);
        breaks.push(None);
    };

    for (start, end, _) in seq_segments.into_iter() {
        let segment = &seq[start..end];
        seqs.push(segment);
        breaks.push(Some(SequenceBreak(end)))
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
    R: TryInto<Builder<3>>,
    I: IntoIterator<Item = Option<R>>,
{
    for (i, (seq, region)) in seq_region_pairs
        .0
        .into_iter()
        .zip(seq_region_pairs.1)
        .enumerate()
    {
        let new_definition = Definition::new(format!("{record_name}_ctg_{i}"), None);
        if let Some(region) = region {
            write_misassembly(
                seq.bytes().collect_vec(),
                std::iter::once(region),
                new_definition,
                writer_fa,
                output_bed.as_mut(),
            )?;
        } else {
            write_misassembly(
                seq.bytes().collect_vec(),
                std::iter::empty::<R>(),
                new_definition,
                writer_fa,
                None,
            )?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_generate_breaks() {
        let seq = "AAAGGCCCGGCCCGGGGATTTTATTTTGGGCCGCCCAATTTAATTT";
        let (seqs, breaks) = generate_breaks(seq, 3, Some(42)).unwrap();
        assert_eq!(
            seqs,
            ["AAAG", "GCCCGGCCCGGGGATTTTATTTTGG", "GCCGCC", "CAATTTAATTT"]
        );
        assert_eq!(
            breaks,
            [
                None,
                Some(SequenceBreak(29)),
                Some(SequenceBreak(35)),
                Some(SequenceBreak(46))
            ]
        );
        assert_eq!(seqs.join(""), seq)
    }
}
