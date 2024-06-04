use std::{collections::HashSet, fs::File, io::Write, ops::Range};

use itertools::Itertools;
use noodles::{
    bed::{self, record::Builder},
    fasta::{
        self,
        record::{Definition, Sequence},
        Writer,
    },
};
use rand::{rngs::StdRng, seq::IteratorRandom, Rng, SeedableRng};

use crate::collapse::Repeat;

/// Generate random sequence segments ranges.
///
/// # Arguments
/// * `seq` - The sequence to generate segments from.
/// * `number` - The number of segments to generate.
/// * `seed` - The random seed to use.
///
/// # Returns
/// An iterator of tuples containing the start, stop, and a random length range starting at the start of the segment.
///
pub fn get_sequence_segments(
    seq: &str,
    number: usize,
    seed: Option<u64>,
) -> Option<impl Iterator<Item = (usize, usize, Range<usize>)>> {
    let mut rng = seed.map_or(StdRng::from_entropy(), StdRng::seed_from_u64);
    let mut positions = (0..seq.len()).choose_multiple(&mut rng, number);
    positions.sort();

    let last_pos = positions.last().copied().map(|p| (p, seq.len()))?;

    Some(
        positions
            .into_iter()
            .tuple_windows()
            .chain(std::iter::once(last_pos))
            .map(move |(start, stop)| {
                let dst = stop - start;
                (start, stop, (start + rng.gen_range(0..dst))..stop)
            }),
    )
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

pub fn flatten_repeats<'a, R: IntoIterator<Item = &'a Repeat>>(seq: &'a str, repeats: R) -> String {
    let mut new_seq = String::with_capacity(seq.len());
    let mut repeats_iter = repeats.into_iter().peekable();

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

pub fn write_misassembly<O, I>(
    seq: Vec<u8>,
    regions: I,
    definition: Definition,
    output_fa: &mut Writer<O>,
    output_bed: Option<&mut bed::Writer<File>>,
) -> eyre::Result<()>
where
    O: Write,
    I: IntoIterator<Item = Builder<3>>,
{
    let full_record_name = definition.to_string();
    let record_name = full_record_name
        .strip_prefix('>')
        .unwrap_or(&full_record_name);

    // Write the BED file if provided.
    if let Some(writer_bed) = output_bed {
        for rp in regions {
            let record = rp.set_reference_sequence_name(record_name).build()?;
            writer_bed.write_record(&record)?;
        }
    };

    output_fa.write_record(&fasta::Record::new(definition, Sequence::from(seq)))?;

    Ok(())
}
