use core::ops::Range;
use iset::{iter::Intervals, IntervalMap};
use itertools::Itertools;
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

pub fn generate_collapse(
    seq: &str,
    repeats: &HashSet<Repeat>,
    num_repeats: usize,
) -> eyre::Result<String> {
    let mut new_seq: String = String::with_capacity(seq.len());

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
    let interval_overlaps: IntervalMap<usize, HashSet<Range<usize>>> =
        IntervalMap::from_iter(intervals.iter(0..seq.len()).map(|(range, repeat)| {
            (
                range.clone(),
                intervals.intervals(range).collect(),
            )
        }));
    println!("Intervals: {:?}", intervals);
    println!("All intervals: {:?}", all_intervals);
    println!("Inteval overlaps: {:?}", interval_overlaps);

    // Iterate thru intervals to construct slices.
    for (range, repeat) in intervals.iter(0..seq.len()) {
        println!("Current: {repeat:?} {range:?}");
        // let (start, end) = (range.start, range.end);

        let complement = interval_overlaps
            .get(range)
            .map(|rs| all_intervals.difference(rs));

        println!("Complement: {:?}", complement);
        // for int in intervals.intervals(range) {
        //     println!("{:?}", int);
        // }

        // // Collapse the repeat.
        // new_seq.push_str(&repeat.seq);

        // // Add the remaining repeats if num_repeats less than repeat count.
        // let remaining_repeats = seq
        //     .get(end - (repeat.seq.len() * repeat.count.saturating_sub(num_repeats))..end)
        //     .unwrap_or_default();

        // new_seq.push_str(remaining_repeats);
    }
    // println!("{:?}", intervals);

    Ok(new_seq)
}

#[cfg(test)]
mod tests {
    use super::*;

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
        let repeats = find_all_repeats(&seq, 5);
        assert_eq!(
            vec![
                Repeat {
                    seq: "TTTTA".to_string(),
                    start: 1,
                    count: 2
                },
                Repeat {
                    seq: "ATTTT".to_string(),
                    start: 0,
                    count: 2
                }
            ],
            repeats.into_iter().collect_vec()
        );
    }

    #[test]
    fn test_find_repeats_multiple() {
        let seq = "GCCCCGCCCCAATTTTAATTTTAATTTT";
        let repeats = find_all_repeats(&seq, 5);
        assert_eq!(
            vec![
                Repeat {
                    seq: "TAATTT".to_string(),
                    start: 15,
                    count: 2
                },
                Repeat {
                    seq: "AATTTT".to_string(),
                    start: 4,
                    count: 3
                },
                Repeat {
                    seq: "AATTTT".to_string(),
                    start: 16,
                    count: 3
                },
                Repeat {
                    seq: "GCCCC".to_string(),
                    start: 0,
                    count: 2
                },
                Repeat {
                    seq: "TTAATT".to_string(),
                    start: 14,
                    count: 2
                },
                Repeat {
                    seq: "TTTTAA".to_string(),
                    start: 12,
                    count: 2
                },
                Repeat {
                    seq: "ATTTTA".to_string(),
                    start: 11,
                    count: 2
                },
                Repeat {
                    seq: "TTTAAT".to_string(),
                    start: 13,
                    count: 2
                }
            ],
            repeats.into_iter().collect_vec()
        )
    }

    #[test]
    fn test_generate_collapse() {
        let seq = "ATTTTATTTT";
        let repeats = find_all_repeats(&seq, 5);
        let new_seq = generate_collapse(seq, &repeats, 20).unwrap();

        assert_eq!("ATTTT", new_seq);
    }

    #[test]
    fn test_generate_collapse_multiple() {
        let seq = "AAAGGCCCGGCCCGGGGATTTTATTTTGGGCCGCCCAATTTAATTT";
        let repeats = find_all_repeats(&seq, 5);
        let new_seq = generate_collapse(seq, &repeats, 20).unwrap();

        // unimplemented!()
    }

    #[test]
    fn test_generate_collapse_interm_seq() {
        let seq = "ATTTTATTTTGCCGAATTTTAATTTTAATTTT";
        let repeats = find_all_repeats(&seq, 5);
        let new_seq = generate_collapse(seq, &repeats, 20).unwrap();
        unimplemented!()
    }
}
