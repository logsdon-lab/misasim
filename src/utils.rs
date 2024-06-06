use std::{fs::File, io::Write, ops::Range};

use eyre::bail;
use iset::{IntervalMap, IntervalSet};
use noodles::{
    bed::{self, record::Builder},
    core::Position,
    fasta::{
        self,
        record::{Definition, Sequence},
        Writer,
    },
};
use rand::{rngs::StdRng, seq::IteratorRandom, SeedableRng};

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
pub fn generate_random_seq_ranges(
    seq_len: usize,
    regions: &IntervalSet<Position>,
    length: usize,
    number: usize,
    seed: Option<u64>,
) -> eyre::Result<Option<impl Iterator<Item = (usize, usize, Range<usize>)>>> {
    let mut rng = seed.map_or(StdRng::from_entropy(), StdRng::seed_from_u64);
    let mut remaining_segments = number;
    let mut positions = IntervalMap::new();

    while remaining_segments > 0 {
        let Some(pos) = regions.unsorted_iter().choose(&mut rng) else {
            break;
        };
        let (start, stop) = (Into::<usize>::into(pos.start), pos.end.into());
        let Some(region_start) = (start..stop).choose(&mut rng) else {
            bail!("Invalid pos: {pos:?}")
        };
        let region_stop = (region_start + 1..region_start + length + 1)
            .choose(&mut rng)
            .map(|stop| stop.clamp(1, seq_len))
            .unwrap();

        // Ensure no overlaps.
        if positions.has_overlap(region_start..region_stop) {
            continue;
        }
        positions.insert(region_start..region_stop, (start, stop));
        remaining_segments -= 1
    }

    let Some(end) = positions.largest().map(|(p, _)| p.end) else {
        bail!("No positions found.")
    };

    // let last_pos = positions.last().context("No positions found.")?;
    Ok(Some(
        positions
            .into_iter(0..end)
            .map(move |(range, (start, stop))| (start, stop, range)),
    ))
}

pub fn write_misassembly<O, R, I>(
    seq: Vec<u8>,
    regions: I,
    definition: Definition,
    output_fa: &mut Writer<O>,
    output_bed: Option<&mut bed::Writer<File>>,
) -> eyre::Result<()>
where
    O: Write,
    R: TryInto<Builder<3>>,
    I: IntoIterator<Item = R>,
{
    let record_name = std::str::from_utf8(definition.name())?;
    // Write the BED file if provided.
    if let Some(writer_bed) = output_bed {
        for builder in regions
            .into_iter()
            .flat_map(|r| TryInto::<Builder<3>>::try_into(r))
        {
            let record = builder.set_reference_sequence_name(record_name).build()?;
            writer_bed.write_record(&record)?;
        }
    };

    output_fa.write_record(&fasta::Record::new(definition, Sequence::from(seq)))?;
    Ok(())
}

#[cfg(test)]
mod test {
    use iset::IntervalSet;
    use itertools::Itertools;
    use noodles::core::Position;

    use super::generate_random_seq_ranges;

    #[test]
    fn test_generate_random_seq_ranges() {
        let positions = vec![Position::new(1).unwrap()..Position::new(10).unwrap()];
        let regions = IntervalSet::from_iter(positions);
        let segments = generate_random_seq_ranges(40, &regions, 10, 2, Some(42))
            .unwrap()
            .unwrap()
            .collect_vec();

        assert_eq!(segments, [(1, 10, 2..3), (1, 10, 3..9)])
    }
}
