use eyre::OptionExt;
use iset::IntervalMap;
use itertools::Itertools;
use rand::prelude::*;

use core::ops::Range;
use std::collections::HashSet;

#[derive(Debug, PartialEq, Eq, Clone, Hash)]
pub struct Repeat {
    pub seq: String,
    pub start: usize,
    pub count: usize,
}

pub fn find_all_repeats(seq: &str, min_len: usize) -> HashSet<Repeat> {
    let rev_seq = seq.chars().rev().collect::<String>();
    let fwd_repeats = find_repeats(seq, min_len);
    let rev_repeats = find_repeats(&rev_seq, min_len);

    // Reorient reverse repeats.
    rev_repeats
        .into_iter()
        .map(|mut repeat| {
            repeat.start = seq.len() - repeat.start - (repeat.seq.len() * repeat.count);
            repeat.seq = repeat.seq.chars().rev().collect();
            repeat
        })
        .chain(fwd_repeats)
        .collect()
}

pub fn find_repeats(seq: &str, min_len: usize) -> Vec<Repeat> {
    let mut repeats: Vec<Repeat> = vec![];
    let mut prev_sfx: Option<(&str, usize)> = None;
    let mut prev_repeat: Option<Repeat> = None;

    for (idx, curr_sfx) in seq
        .char_indices()
        .flat_map(|(i, _)| seq.get(i..).map(|s| (i, s)))
        .sorted_by(|a, b| a.1.cmp(b.1))
    {
        if let Some((mut sfx, sfx_pos)) = prev_sfx {
            // Does current suffix start with this suffix?
            if !curr_sfx.starts_with(sfx) {
                // Remove prev suffix from current suffix.
                // This isolates repeated elements at start of current suffix.
                prev_sfx = curr_sfx
                    .strip_suffix(sfx)
                    .map_or(Some((curr_sfx, 0)), |s| Some((s, 0)));

                if let Some((s, _)) = prev_sfx {
                    sfx = s;
                }
                if let Some(repeat) = prev_repeat {
                    repeats.push(repeat);
                    prev_repeat = None;
                }
            }
            let prev_sfx_len = sfx.len();
            // Remove small suffixes.
            if prev_sfx_len < min_len {
                prev_sfx = None;
                continue;
            }
            // Check slice ahead.
            let sfx_end = sfx_pos + prev_sfx_len;
            let Some(next_sfx) = &curr_sfx.get(sfx_end..sfx_end + prev_sfx_len) else {
                continue;
            };

            // If same as previous suffix, increment count of repeat.
            // Otherwise, store repeat.
            if *next_sfx == sfx {
                if let Some(repeat) = prev_repeat.as_mut() {
                    repeat.count += 1;
                } else {
                    prev_repeat = Some(Repeat {
                        seq: sfx.to_owned(),
                        start: idx,
                        count: 2,
                    });
                }
            } else {
                continue;
            }
            prev_sfx = Some((sfx, sfx_end));
        } else {
            prev_sfx = Some((curr_sfx, 0));
        }
    }

    if let Some(repeat) = prev_repeat {
        repeats.push(repeat);
    }

    repeats
}

pub fn flatten_repeats<'a, R: IntoIterator<Item = &'a Repeat>>(
    seq: &str,
    repeats: R,
    rng: &mut StdRng,
    num_repeats: usize,
) -> String {
    let mut new_seq = String::with_capacity(seq.len());
    // Choose a random number of repeats to collapse
    let repeats = repeats
        .into_iter()
        .sorted_by(|a, b| a.start.cmp(&b.start))
        .choose_multiple(rng, num_repeats);
    let mut repeats_iter = repeats.iter().peekable();

    // Get first segment before repeat.
    if let Some(first_segment) = repeats_iter
        .peek()
        .map(|r| seq.get(..r.start))
        .unwrap_or_default()
    {
        new_seq.push_str(first_segment);
    }

    while let Some(rp) = repeats_iter.next() {
        let rp_end = rp.start + (rp.seq.len() * rp.count);
        // Collapse repeat.
        new_seq.push_str(&rp.seq);

        // Add the segment between this repeat and the next one.
        let next_rp = repeats_iter.peek();
        if let Some(next_repeat) = next_rp.and_then(|r| seq.get(rp_end..r.start)) {
            new_seq.push_str(next_repeat);
        } else if let Some(last_segment) = seq.get(rp_end..) {
            new_seq.push_str(last_segment);
        }
    }

    new_seq
}

pub fn generate_collapse(
    seq: &str,
    repeats: &HashSet<Repeat>,
    num_repeats: usize,
    seed: Option<u64>,
) -> eyre::Result<String> {
    let mut rng = seed.map_or(StdRng::from_entropy(), StdRng::seed_from_u64);
    let mut new_seqs: Vec<String> = vec![];
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

        let new_seq = flatten_repeats(
            seq,
            complement
                .into_iter()
                .flat_map(|r| intervals.get(r.clone()))
                .chain(std::iter::once(repeat)),
            &mut rng,
            num_repeats,
        );

        new_seqs.push(new_seq);
    }

    // If no non-overlapping intervals, use only unique intervals.
    if new_seqs.is_empty() {
        let new_seq = flatten_repeats(
            seq,
            intervals.iter(0..seq.len()).map(|(_, r)| r),
            &mut rng,
            num_repeats,
        );
        new_seqs.push(new_seq);
    }

    // Choose a random new sequence.
    new_seqs
        .choose(&mut rng)
        .ok_or_eyre("No collapsed sequences can be generated.")
        .cloned()
}

#[cfg(test)]
mod tests {
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
        let repeats = find_all_repeats(&seq, 5);
        let new_seq = generate_collapse(seq, &repeats, 20, None).unwrap();

        assert_eq!("ATTTT", new_seq);
    }

    #[test]
    fn test_generate_collapse_multiple() {
        let seq = "AAAGGCCCGGCCCGGGGATTTTATTTTGGGCCGCCCAATTTAATTT";
        let repeats = find_all_repeats(&seq, 5);
        let new_seq = generate_collapse(seq, &repeats, 4, Some(42)).unwrap();

        assert_eq!(new_seq, "AAAGGCCCGGGGATTTTGGGCCGCCCAATTT")
    }
}
