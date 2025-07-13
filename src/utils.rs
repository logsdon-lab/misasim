use eyre::bail;
use itertools::Itertools;
use rand::{rngs::StdRng, seq::IteratorRandom, SeedableRng};
use rust_lapper::{Interval, Lapper};
use std::{fmt::Debug, ops::Range};

use crate::sequence::{SequenceSegment, SequenceType};

/// Generate random sequence segments ranges.
///
/// # Arguments
/// * `seq_len` - Length of sequence to clamp regions.
/// * `regions` - Positions to choose segments from.
/// * `length` - The maximum length of a generate segment.
/// * `number` - The number of segments to generate.
/// * `seed` - The random seed to use.
///
/// # Returns
/// An iterator of tuples containing the start, stop, and a random length range starting at the start of the segment.
///
pub fn generate_random_seq_ranges<T>(
    seq_len: usize,
    regions: &Lapper<usize, T>,
    length: usize,
    number: usize,
    seed: Option<u64>,
    randomize_length: bool,
) -> eyre::Result<impl Iterator<Item = Interval<usize, Range<usize>>>>
where
    T: Eq + Clone + Send + Sync + Debug,
{
    let mut rng = seed.map_or(StdRng::from_entropy(), StdRng::seed_from_u64);
    let mut remaining_segments = number;
    let mut positions = Lapper::new(vec![]);
    // Keep going until required number of segments generated
    while remaining_segments > 0 {
        // Choose a starting position within the provided region set. ex. bed file.
        let Some(pos) = regions.iter().choose(&mut rng) else {
            break;
        };
        let (start, stop) = (pos.start, pos.stop);
        // Then if randomizing length, choose a starting position within the selected region.
        // Choose a random ending position.
        let (region_start, region_stop) = if randomize_length {
            let Some(region_start) = (start..stop).choose(&mut rng) else {
                bail!("Invalid pos: {pos:?}")
            };
            let region_stop = (region_start + 1..region_start + length + 1)
                .choose(&mut rng)
                .map(|stop| stop.clamp(1, seq_len))
                .unwrap();
            (region_start, region_stop)
        } else {
            // Choose a starting position within the range shortened by the desired length.
            // Use the randomly selected starting position and add the length.
            let stop = stop - length;
            let Some(region_start) = (start..stop).choose(&mut rng) else {
                bail!("Invalid pos: {pos:?}")
            };
            (region_start, region_start + length)
        };

        // Ensure no overlaps.
        // Keep iterating until a valid position found.
        if positions.find(region_start, region_stop).count() != 0 {
            continue;
        }
        positions.insert(Interval {
            start: region_start,
            stop: region_stop,
            val: start..stop,
        });
        remaining_segments -= 1
    }

    // let last_pos = positions.last().context("No positions found.")?;
    Ok(positions.into_iter().map(|itv| Interval {
        start: itv.val.start,
        stop: itv.val.end,
        val: itv.start..itv.stop,
    }))
}

/// Subtract interval by a list of non-overlapping intervals.
pub fn subtract_intervals<T>(
    itv: Interval<usize, T>,
    other: impl Iterator<Item = Interval<usize, T>>,
) -> Vec<Interval<usize, T>>
where
    T: Eq + Clone + Send + Sync,
{
    let mut split_intervals = Vec::new();
    let mut st = itv.start;
    let mut last = itv.stop;
    for ovl_itv in other.into_iter().sorted_by(|a, b| a.start.cmp(&b.start)) {
        if last >= ovl_itv.start && last <= ovl_itv.stop {
            //    |---|
            // * |---|
            last = ovl_itv.start;
        } else if st <= ovl_itv.stop && st >= ovl_itv.start {
            //   |---|
            // *  |---|
            st = ovl_itv.stop;
        } else if st >= ovl_itv.start && last <= ovl_itv.stop {
            //   |---|
            // * |---|
            break;
        } else if ovl_itv.start > st && ovl_itv.stop < last {
            //    |-|
            // * |---|
            split_intervals.push(Interval {
                start: st,
                stop: ovl_itv.start,
                val: itv.val.clone(),
            });
            st = ovl_itv.stop;
        }
    }
    // Add remainder.
    if st != last {
        split_intervals.push(Interval {
            start: st,
            stop: last,
            val: itv.val,
        });
    }
    split_intervals
}

/// Subtract misassemblies from a sequence.
/// # Args
/// * `seq`
///     * Complete sequence.
///     * `[0, seq.len())`
/// * `misassemblies`
///     * `SequenceSegments` iterator within `seq` coordinates.
///
/// # Returns:
/// * Good intervals with `None` val.
pub fn subtract_misassembled_sequences<'a>(
    seq: &str,
    misassemblies: impl Iterator<Item = &'a SequenceSegment>,
) -> Vec<SequenceSegment> {
    let mut split_intervals = Vec::new();
    let mut st = 0;
    let mut last = seq.len();
    for misassembly in misassemblies
        .into_iter()
        .sorted_by(|a, b| a.itv.start.cmp(&b.itv.start))
    {
        if last >= misassembly.itv.start && last <= misassembly.itv.stop {
            //    |---|
            // * |---|
            last = misassembly.itv.start;
        } else if st <= misassembly.itv.stop && st >= misassembly.itv.start {
            //   |---|
            // *  |---|
            st = misassembly.itv.stop;
        } else if st >= misassembly.itv.start && last <= misassembly.itv.stop {
            //   |---|
            // * |---|
            break;
        } else if misassembly.itv.start > st && misassembly.itv.stop < last {
            //    |-|
            // * |---|
            split_intervals.push(SequenceSegment {
                typ: SequenceType::Good,
                itv: Interval {
                    start: st,
                    stop: misassembly.itv.start,
                    val: None,
                },
            });
            st = misassembly.itv.stop;
        }
    }
    // Add remainder.
    if st != last {
        split_intervals.push(SequenceSegment {
            typ: SequenceType::Good,
            itv: Interval {
                start: st,
                stop: last,
                val: None,
            },
        });
    }
    split_intervals
}

pub fn calculate_new_coords(seqs: &[SequenceSegment]) -> Vec<Interval<usize, ()>> {
    let mut adj_coords = Vec::with_capacity(seqs.len());
    let mut delta: isize = 0;
    for seq in seqs {
        // Adjust coordinates of coordinates.
        let (new_start, new_stop) = if delta.is_negative() {
            let delta_usize = -delta as usize;
            (seq.itv.start - delta_usize, seq.itv.stop - delta_usize)
        } else {
            let delta_usize = delta as usize;
            (seq.itv.start + delta_usize, seq.itv.stop + delta_usize)
        };
        match seq.typ {
            SequenceType::Good
            | SequenceType::Gap
            | SequenceType::Break
            | SequenceType::Inversion => {
                adj_coords.push(Interval {
                    start: new_start,
                    stop: new_stop,
                    val: (),
                });
            }
            SequenceType::Misjoin => {
                let adj_delta = seq.itv.stop - seq.itv.start;
                delta -= adj_delta as isize;
                // Deleted from assembly. Null interval
                adj_coords.push(Interval {
                    start: new_start,
                    stop: new_start,
                    val: (),
                });
            }
            SequenceType::FalseDuplication => {
                let dupe_seq = seq
                    .itv
                    .val
                    .as_ref()
                    .expect("Invalid state. False dupe with no sequence.");
                let adj_delta = dupe_seq.len() - (seq.itv.stop - seq.itv.start);
                delta += adj_delta as isize;
                // Add duplicate sequence length to end to match.
                adj_coords.push(Interval {
                    start: new_start,
                    stop: new_stop + adj_delta,
                    val: (),
                });
            }
        }
    }

    adj_coords
}

#[cfg(test)]
mod test {
    use itertools::Itertools;
    use rust_lapper::{Interval, Lapper};

    use crate::{
        sequence::{SequenceSegment, SequenceType},
        utils::calculate_new_coords,
    };

    use super::generate_random_seq_ranges;

    #[test]
    fn test_generate_random_seq_ranges() {
        let regions = Lapper::new(vec![Interval {
            start: 1,
            stop: 10,
            val: (),
        }]);
        let segments = generate_random_seq_ranges(40, &regions, 10, 2, Some(42), true)
            .unwrap()
            .collect_vec();

        assert_eq!(
            segments,
            [
                Interval {
                    start: 1,
                    stop: 10,
                    val: 2..3
                },
                Interval {
                    start: 1,
                    stop: 10,
                    val: 3..9
                }
            ]
        )
    }

    #[test]
    fn test_generate_random_seq_ranges_static_length() {
        let regions = Lapper::new(vec![Interval {
            start: 1,
            stop: 10,
            val: (),
        }]);
        // Generate two regions of length 2.
        let segments = generate_random_seq_ranges(40, &regions, 2, 2, Some(42), false)
            .unwrap()
            .collect_vec();
        assert_eq!(
            segments,
            [
                Interval {
                    start: 1,
                    stop: 10,
                    val: 4..6
                },
                Interval {
                    start: 1,
                    stop: 10,
                    val: 7..9
                }
            ]
        )
    }

    #[test]
    fn subtract_misassembled_sequences() {}

    #[test]
    fn test_compute_delta_seq_segments() {
        // let seq = "ATTATTATTGCA";
        let seqs = vec![
            // ATT---------
            SequenceSegment {
                typ: SequenceType::Good,
                itv: Interval {
                    start: 0,
                    stop: 3,
                    val: None,
                },
            },
            // ---AT-------
            SequenceSegment {
                typ: SequenceType::FalseDuplication,
                itv: Interval {
                    start: 3,
                    stop: 5,
                    val: Some("ATAT".to_owned()),
                },
            },
            // -----TAT----
            SequenceSegment {
                typ: SequenceType::Misjoin,
                itv: Interval {
                    start: 5,
                    stop: 7,
                    val: None,
                },
            },
            // --------TGCA
            SequenceSegment {
                typ: SequenceType::Good,
                itv: Interval {
                    start: 7,
                    stop: 12,
                    val: None,
                },
            },
        ];
        let new_coords = calculate_new_coords(&seqs);
        assert_eq!(
            new_coords,
            vec![
                Interval {
                    start: 0,
                    stop: 3,
                    val: ()
                },
                // This gets 2 additional bases from duplication.
                Interval {
                    start: 3,
                    stop: 7,
                    val: ()
                },
                // This is gone.
                Interval {
                    start: 7,
                    stop: 7,
                    val: ()
                },
                Interval {
                    start: 7,
                    stop: 12,
                    val: ()
                },
            ]
        )
    }
}
