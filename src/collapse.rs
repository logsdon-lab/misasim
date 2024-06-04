use core::ops::Range;
use std::collections::HashSet;

use eyre::OptionExt;
use iset::IntervalMap;
use itertools::Itertools;
use noodles::{
    bed::{
        self,
        record::{Builder, OptionalFields},
    },
    core::Position,
};
use rand::prelude::*;

use crate::utils::{find_all_repeats, flatten_repeats};

#[derive(Debug, PartialEq, Eq, Clone, Hash)]
pub struct Repeat {
    pub seq: String,
    pub start: usize,
    pub count: usize,
}

impl From<Repeat> for Builder<3> {
    fn from(rp: Repeat) -> Self {
        bed::Record::<3>::builder()
            .set_start_position(Position::new(rp.start.clamp(1, usize::MAX)).unwrap())
            .set_end_position(Position::new(rp.start + (rp.seq.len() * rp.count)).unwrap())
            .set_optional_fields(OptionalFields::from(vec![rp.count.to_string(), rp.seq]))
    }
}

/// Collapsed sequence and their repeats
// TODO: Add function to reform full sequence.
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct CollapsedSequence {
    pub seq: String,
    pub repeats: Vec<Repeat>,
}

pub fn generate_collapse(
    seq: &str,
    min_length: usize,
    num_repeats: usize,
    seed: Option<u64>,
) -> eyre::Result<CollapsedSequence> {
    let mut rng = seed.map_or(StdRng::from_entropy(), StdRng::seed_from_u64);
    let mut new_seqs: Vec<CollapsedSequence> = vec![];
    let repeats = find_all_repeats(seq, min_length);
    let intervals: IntervalMap<usize, Repeat> =
        IntervalMap::from_iter(repeats.iter().map(|repeat| {
            let (start, stop) = (
                repeat.start,
                repeat.start + (repeat.seq.len() * repeat.count),
            );
            (start..stop, repeat.clone())
        }));
    let all_intervals: HashSet<Range<usize>> = intervals
        .iter(0..seq.len())
        .map(|(range, _)| range)
        .collect();
    let interval_overlaps: IntervalMap<usize, HashSet<Range<usize>>> = IntervalMap::from_iter(
        intervals
            .iter(0..seq.len())
            .map(|(range, _)| (range.clone(), intervals.intervals(range).collect())),
    );

    // Iterate thru intervals to construct slices.
    for (range, repeat) in intervals
        .iter(0..seq.len())
        // Skip unique, non-overlapping intervals.
        .filter(|(r, _)| {
            interval_overlaps
                .get(r.clone())
                .map(|i| i.len() != 1)
                .unwrap_or_default()
        })
    {
        let Some(complement) = interval_overlaps
            .get(range)
            .map(|rs| all_intervals.difference(rs))
        else {
            continue;
        };

        let final_repeats = complement
            .into_iter()
            .flat_map(|r| intervals.get(r.clone()))
            .chain(std::iter::once(repeat))
            .cloned()
            .sorted_by(|a, b| a.start.cmp(&b.start))
            .choose_multiple(&mut rng, num_repeats);
        let new_seq = flatten_repeats(seq, final_repeats.iter());

        new_seqs.push(CollapsedSequence {
            seq: new_seq,
            repeats: final_repeats,
        });
    }

    // If no non-overlapping intervals, use only unique intervals.
    if new_seqs.is_empty() {
        let final_repeats = intervals
            .iter(0..seq.len())
            .map(|(_, r)| r)
            .cloned()
            .sorted_by(|a, b| a.start.cmp(&b.start))
            .choose_multiple(&mut rng, num_repeats);
        let new_seq = flatten_repeats(seq, final_repeats.iter());
        new_seqs.push(CollapsedSequence {
            seq: new_seq,
            repeats: final_repeats,
        });
    }

    // Choose a random new sequence.
    new_seqs
        .into_iter()
        .choose(&mut rng)
        .ok_or_eyre("No collapsed sequences can be generated.")
}

#[cfg(test)]
mod tests {
    use crate::utils::find_all_repeats;

    use super::*;

    fn sort_repeats(repeats: &mut Vec<Repeat>) {
        repeats.sort_by(|a: &Repeat, b: &Repeat| a.start.cmp(&b.start));
    }

    #[test]
    fn test_find_repeats() {
        let seq = "ATTTTATTTT";
        let repeats = find_all_repeats(&seq, 5);
        assert_eq!(
            vec![Repeat {
                seq: String::from("ATTTT"),
                start: 0,
                count: 2
            }],
            repeats.into_iter().collect_vec()
        );
    }

    #[test]
    fn test_find_repeats_overlap() {
        let seq = "ATTTTATTTTA";
        let mut repeats = find_all_repeats(&seq, 5).into_iter().collect_vec();
        let exp_repeats = vec![
            Repeat {
                seq: "ATTTT".to_string(),
                start: 0,
                count: 2,
            },
            Repeat {
                seq: "TTTTA".to_string(),
                start: 1,
                count: 2,
            },
        ];
        sort_repeats(&mut repeats);
        assert_eq!(exp_repeats, repeats);
    }

    #[test]
    fn test_find_repeats_multiple() {
        let seq = "GCCCCGCCCCAATTTTAATTTTAATTTT";
        let mut repeats = find_all_repeats(&seq, 5).into_iter().collect_vec();
        let mut exp_repeats = vec![
            Repeat {
                seq: "TAATTT".to_string(),
                start: 15,
                count: 2,
            },
            Repeat {
                seq: "AATTTT".to_string(),
                start: 4,
                count: 3,
            },
            Repeat {
                seq: "AATTTT".to_string(),
                start: 16,
                count: 3,
            },
            Repeat {
                seq: "GCCCC".to_string(),
                start: 0,
                count: 2,
            },
            Repeat {
                seq: "TTAATT".to_string(),
                start: 14,
                count: 2,
            },
            Repeat {
                seq: "TTTTAA".to_string(),
                start: 12,
                count: 2,
            },
            Repeat {
                seq: "ATTTTA".to_string(),
                start: 11,
                count: 2,
            },
            Repeat {
                seq: "TTTAAT".to_string(),
                start: 13,
                count: 2,
            },
        ];
        sort_repeats(&mut repeats);
        sort_repeats(&mut exp_repeats);
        assert_eq!(exp_repeats, exp_repeats)
    }

    #[test]
    fn test_generate_collapse() {
        let seq = "ATTTTATTTT";
        let new_seq = generate_collapse(seq, 5, 20, None).unwrap();
        assert_eq!(
            CollapsedSequence {
                seq: "ATTTT".to_string(),
                repeats: [Repeat {
                    seq: "ATTTT".to_string(),
                    start: 0,
                    count: 2
                }]
                .to_vec()
            },
            new_seq
        );
    }

    #[test]
    fn test_generate_collapse_multiple() {
        let seq = "AAAGGCCCGGCCCGGGGATTTTATTTTGGGCCGCCCAATTTAATTT";
        let new_seq = generate_collapse(seq, 5, 4, Some(42)).unwrap();

        assert_eq!(
            CollapsedSequence {
                seq: "AAAGGCCCGGGGATTTTGGGCCGCCCAATTT".to_string(),
                repeats: [
                    Repeat {
                        seq: "GGCCC".to_string(),
                        start: 3,
                        count: 2
                    },
                    Repeat {
                        seq: "ATTTT".to_string(),
                        start: 17,
                        count: 2
                    },
                    Repeat {
                        seq: "AATTT".to_string(),
                        start: 36,
                        count: 2
                    }
                ]
                .to_vec()
            },
            new_seq
        );
    }
}
