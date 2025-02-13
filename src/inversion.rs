use eyre::ContextCompat;
use iset::IntervalSet;
use itertools::Itertools;
use noodles::{
    bed::{
        record::{Builder, OptionalFields},
        Record,
    },
    core::Position,
};

use crate::utils::generate_random_seq_ranges;

#[derive(Debug, PartialEq, Eq, Clone)]
pub struct RevCompSequence {
    pub start: usize,
    pub end: usize,
    pub seq: String,
}

impl TryFrom<RevCompSequence> for Builder<3> {
    type Error = eyre::Error;

    fn try_from(rev_seq: RevCompSequence) -> Result<Self, eyre::Error> {
        Ok(Record::builder()
            .set_start_position(Position::new(rev_seq.start).context("Zero start position")?)
            .set_end_position(Position::new(rev_seq.end).context("Zero end position")?)
            .set_optional_fields(OptionalFields::from(vec![rev_seq.seq])))
    }
}

#[derive(Debug, PartialEq, Eq)]
pub struct InvertedSequence {
    pub seq: String,
    pub inverted_seqs: Vec<RevCompSequence>,
}

pub fn generate_inversion(
    seq: &str,
    regions: &IntervalSet<Position>,
    length: usize,
    number_inv: usize,
    seed: Option<u64>,
    randomize_length: bool,
) -> eyre::Result<InvertedSequence> {
    let mut new_seq = String::with_capacity(seq.len());
    let mut inverted_seqs: Vec<RevCompSequence> = Vec::with_capacity(number_inv);
    let seq_segments = generate_random_seq_ranges(
        seq.len(),
        regions,
        length,
        number_inv,
        seed,
        randomize_length,
    )?
    .context("No sequence segments")?
    .collect_vec();

    let mut seq_iter = seq_segments.into_iter().peekable();
    // Add starting sequence before first position.
    if let Some((_, _, inv_range)) = seq_iter.peek() {
        new_seq.push_str(&seq[..inv_range.start]);
    };

    while let Some((_, _, rrange)) = seq_iter.next() {
        let inv_seq: String = seq[rrange.clone()]
            .chars()
            .map(|nt| match nt.to_ascii_uppercase() {
                'A' => 'T',
                'T' => 'A',
                'G' => 'C',
                'C' => 'G',
                _ => panic!("Invalid {nt}"),
            })
            .rev()
            .collect();

        new_seq.push_str(&inv_seq);

        inverted_seqs.push(RevCompSequence {
            start: rrange.start,
            end: rrange.end,
            seq: inv_seq,
        });

        let remaining_seq = if let Some((_, _, next_rrange)) = seq_iter.peek() {
            &seq[rrange.end..next_rrange.start]
        } else {
            &seq[rrange.end..seq.len()]
        };

        new_seq.push_str(remaining_seq);
    }

    Ok(InvertedSequence {
        seq: new_seq,
        inverted_seqs,
    })
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_generate_inversion() {
        let seq = "AAAGGCCCGGCCCGGGGATTTTATTTTGGGCCGCCCAATTTAATTT";
        let regions = IntervalSet::from_iter(std::iter::once(
            Position::new(1).unwrap()..Position::new(seq.len()).unwrap(),
        ));
        let new_seq = generate_inversion(seq, &regions, 10, 1, Some(100), true).unwrap();

        assert_eq!(
            InvertedSequence {
                seq: "AAAGGCCCGGCCCGGGGATTTTATTTTGGGCCGCCCAAATTAAATT".to_string(),
                inverted_seqs: [RevCompSequence {
                    start: 38,
                    end: 44,
                    seq: "ATTAAA".to_string()
                }]
                .to_vec()
            },
            new_seq
        );
    }
}
