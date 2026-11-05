use crate::{format::SeedMask, io::Int};
use anyhow::{anyhow, Result};
use log::info;
use rayon::prelude::*;
use std::{ops::Range, sync::Mutex, time::Instant};

/// The most bucket numbers a text is partitioned into.
const MAX_BUCKETS: usize = 1 << 16;

/// The suffixes one task of a pass over the text takes.
const CHUNK_SIZE: usize = 1 << 19;

/// The rank of each byte of the text in byte order, and how many leading
/// characters a bucket number is made of.
pub(super) struct PrefixEncoder {
    table: [u8; 256],

    /// The number of ranks, which is the base a bucket number is
    /// written in.
    radix: usize,

    /// The number of leading characters in a bucket number.
    count: usize,

    /// The weight of the top digit, `radix^(count-1)`.
    top: usize,

    buckets: usize,
}

impl PrefixEncoder {
    /// The encoder for a text whose bytes `occupancy` marks.
    pub(super) fn from_occupancy(occupancy: [bool; 256]) -> Self {
        let mut table = [0u8; 256];
        let mut radix = 0;
        for (c, &seen) in occupancy.iter().enumerate() {
            if seen {
                table[c] = radix as u8;
                radix += 1;
            }
        }

        // as many characters as keep the bucket numbers under
        // MAX_BUCKETS; the cap is for a radix of 1, which never
        // reaches it
        let mut count = 0;
        let mut buckets = 1;
        while count < 16 && buckets * radix <= MAX_BUCKETS {
            buckets *= radix;
            count += 1;
        }
        info!("{radix} distinct bytes: {count}-character prefix, {buckets} buckets");

        Self {
            table,
            radix,
            count,
            top: buckets / radix.max(1),
            buckets,
        }
    }

    /// The number of distinct bucket numbers.
    pub(super) fn num_buckets(&self) -> usize {
        self.buckets
    }

    /// The bucket number of the suffix at `i`: the ranks of its first
    /// `count` characters, or of those at the care positions of `mask`,
    /// as the digits of a base-`radix` number.
    #[inline(always)]
    pub(super) fn encode(&self, text: &[u8], i: usize, mask: Option<&SeedMask>) -> u16 {
        let mut key = 0;
        for j in 0..self.count {
            // a character past the end of the text or the mask
            // ranks 0, so a short suffix keeps every digit
            let c = match mask {
                Some(mask) => mask.positions.get(j).and_then(|&off| text.get(i + off)),
                None => text.get(i + j),
            };
            let rank = c.map_or(0, |&c| self.table[c as usize] as usize);
            key = key * self.radix + rank;
        }
        key as u16
    }

    /// The bucket number of the suffix after the one at `i` with
    /// bucket number `key`, without a mask.
    #[inline(always)]
    pub(super) fn roll(&self, text: &[u8], i: usize, key: usize) -> usize {
        // the next suffix drops this one's first character from the
        // top digit and takes the character `count` ahead as its last:
        //
        //     key(i+1) = (key(i) - rank(text[i]) × radix^(count-1)) × radix
        //              + rank(text[i+count])
        let out = self.table[text[i] as usize] as usize;
        let incoming = text
            .get(i + self.count)
            .map_or(0, |&c| self.table[c as usize] as usize);
        (key - out * self.top) * self.radix + incoming
    }
}

/// The partitions: runs of whole buckets in prefix order.
pub(super) struct Layout {
    /// The partition of each bucket.
    partition_of: Vec<usize>,

    /// The buckets of each partition.
    pub(super) buckets: Vec<Range<usize>>,

    /// The suffixes of each partition, as a range of the suffix array.
    pub(super) ranges: Vec<Range<usize>>,
}

impl Layout {
    /// Count the suffixes of `text` in each bucket, those whose first
    /// byte `keep` marks, and cut the buckets into at most
    /// `num_partitions` partitions of about the same number of
    /// suffixes.
    pub(super) fn new(
        text: &[u8],
        encoder: &PrefixEncoder,
        keep: &[bool; 256],
        mask: Option<&SeedMask>,
        num_partitions: usize,
    ) -> Self {
        let now = Instant::now();
        let num_buckets = encoder.num_buckets();
        let counts = (0..text.len().div_ceil(CHUNK_SIZE))
            .into_par_iter()
            .fold(
                || vec![0; num_buckets],
                |mut counts, chunk| {
                    let start = chunk * CHUNK_SIZE;
                    let end = (start + CHUNK_SIZE).min(text.len());
                    match mask {
                        Some(mask) => {
                            for i in start..end {
                                if keep[text[i] as usize] {
                                    let key = encoder.encode(text, i, Some(mask));
                                    counts[key as usize] += 1;
                                }
                            }
                        }
                        None => {
                            let mut key = encoder.encode(text, start, None) as usize;
                            for i in start..end {
                                if keep[text[i] as usize] {
                                    counts[key] += 1;
                                }
                                key = encoder.roll(text, i, key);
                            }
                        }
                    }
                    counts
                },
            )
            .reduce(
                || vec![0; num_buckets],
                |mut total, counts| {
                    for (total, count) in total.iter_mut().zip(counts) {
                        *total += count;
                    }
                    total
                },
            );

        let layout = Self::cut(&counts, num_partitions);
        info!(
            "Counted {} suffixes into {num_buckets} buckets and {} partitions in {:?}",
            layout.num_suffixes(),
            layout.ranges.len(),
            now.elapsed()
        );
        layout
    }

    fn cut(counts: &[usize], num_partitions: usize) -> Self {
        let num_suffixes: usize = counts.iter().sum();
        let share = num_suffixes.div_ceil(num_partitions.max(1));
        let mut partition_of = vec![0; counts.len()];
        let mut buckets = vec![];
        let mut ranges: Vec<Range<usize>> = vec![];
        let mut first_bucket = 0;
        let mut start = 0;
        let mut taken = 0;
        for (bucket, &count) in counts.iter().enumerate() {
            partition_of[bucket] = ranges.len();
            taken += count;

            // a partition closes once the suffixes taken so far
            // reach its cumulative share, so a bucket far over
            // its share is paid back by the partitions after it
            if taken > start && taken >= share * (ranges.len() + 1) {
                buckets.push(first_bucket..bucket + 1);
                ranges.push(start..taken);
                first_bucket = bucket + 1;
                start = taken;
            }
        }
        if taken > start {
            buckets.push(first_bucket..counts.len());
            ranges.push(start..taken);
        }
        Self {
            partition_of,
            buckets,
            ranges,
        }
    }

    pub(super) fn num_suffixes(&self) -> usize {
        self.ranges.last().map_or(0, |range| range.end)
    }

    /// Runs of partitions of at most `max_suffixes` suffixes each, or
    /// one partition where that alone is more.
    pub(super) fn turns(&self, max_suffixes: usize) -> Vec<Range<usize>> {
        let mut turns = vec![];
        let mut first = 0;
        let mut taken = 0;
        for (p, range) in self.ranges.iter().enumerate() {
            if taken > 0 && taken + range.len() > max_suffixes {
                turns.push(first..p);
                first = p;
                taken = 0;
            }
            taken += range.len();
        }
        if first < self.ranges.len() {
            turns.push(first..self.ranges.len());
        }
        turns
    }
}

/// `values` in the order of `key(i)` for each index `i`, every key
/// below `num_keys`, and where each key's run starts, with one more
/// entry for the end.
pub(super) fn counting_sort<T: Copy + Default>(
    values: &[T],
    num_keys: usize,
    key: impl Fn(usize) -> usize,
) -> (Vec<T>, Vec<usize>) {
    let mut starts = vec![0; num_keys + 1];
    for i in 0..values.len() {
        starts[key(i) + 1] += 1;
    }
    for k in 0..num_keys {
        starts[k + 1] += starts[k];
    }
    let mut next = starts.clone();
    let mut ordered = vec![T::default(); values.len()];
    for (i, &value) in values.iter().enumerate() {
        let slot = &mut next[key(i)];
        ordered[*slot] = value;
        *slot += 1;
    }
    (ordered, starts)
}

/// The suffixes of each partition of `partitions` with their buckets,
/// from one scan of `text` keeping those whose first byte `keep` marks.
pub(super) fn gather<W: Int>(
    text: &[u8],
    encoder: &PrefixEncoder,
    keep: &[bool; 256],
    mask: Option<&SeedMask>,
    layout: &Layout,
    partitions: Range<usize>,
) -> Result<Vec<(Vec<W>, Vec<u16>)>> {
    let now = Instant::now();
    let first = partitions.start;
    let count = partitions.len();
    // one compare for "in this turn's buckets", so the branch is
    // rarely taken whichever turn this is; as lo <= key && key < hi
    // the first test is a coin flip in the middle turns and the
    // scan ran three times slower there
    let lo = layout.buckets[first].start;
    let width = layout.buckets[partitions.end - 1].end - lo;
    let gathered: Vec<Mutex<(Vec<W>, Vec<u16>)>> = partitions
        .map(|p| {
            let len = layout.ranges[p].len();
            Mutex::new((Vec::with_capacity(len), Vec::with_capacity(len)))
        })
        .collect();

    // the buffers are sized once per thread and filled by index,
    // so the loop holds no call and keeps its state in registers
    (0..text.len().div_ceil(CHUNK_SIZE))
        .into_par_iter()
        .try_for_each_init(
            || (vec![W::default(); CHUNK_SIZE], vec![0u16; CHUNK_SIZE]),
            |(suffixes, buckets), chunk| -> Result<()> {
                let start = chunk * CHUNK_SIZE;
                let end = (start + CHUNK_SIZE).min(text.len());
                let mut n = 0;
                match mask {
                    Some(mask) => {
                        for i in start..end {
                            let key = encoder.encode(text, i, Some(mask)) as usize;
                            if keep[text[i] as usize] && key.wrapping_sub(lo) < width {
                                suffixes[n] = W::from_usize(i);
                                buckets[n] = key as u16;
                                n += 1;
                            }
                        }
                    }
                    None => {
                        let mut key = encoder.encode(text, start, None) as usize;
                        for i in start..end {
                            if keep[text[i] as usize] && key.wrapping_sub(lo) < width {
                                suffixes[n] = W::from_usize(i);
                                buckets[n] = key as u16;
                                n += 1;
                            }
                            key = encoder.roll(text, i, key);
                        }
                    }
                }
                let suffixes = &suffixes[..n];
                let buckets = &buckets[..n];

                let partition =
                    |i: usize| layout.partition_of[buckets[i] as usize] - first;
                let (suffixes, starts) = counting_sort(suffixes, count, partition);
                let (buckets, _) = counting_sort(buckets, count, partition);
                for (p, run) in starts.windows(2).enumerate() {
                    let run = run[0]..run[1];
                    if !run.is_empty() {
                        let mut g = gathered[p].lock().map_err(|e| anyhow!("{e}"))?;
                        g.0.extend_from_slice(&suffixes[run.clone()]);
                        g.1.extend_from_slice(&buckets[run]);
                    }
                }
                Ok(())
            },
        )?;

    let gathered = gathered
        .into_iter()
        .map(|g| g.into_inner().map_err(|e| anyhow!("{e}")))
        .collect::<Result<Vec<_>>>()?;
    info!(
        "Gathered {} suffixes of partitions {first}..{} in {:?}",
        gathered.iter().map(|g| g.0.len()).sum::<usize>(),
        first + count,
        now.elapsed()
    );
    Ok(gathered)
}

#[cfg(test)]
mod test {
    use super::{counting_sort, gather, Layout, PrefixEncoder};
    use crate::format::SeedMask;
    use anyhow::Result;
    use pretty_assertions::assert_eq;

    fn encoder(text: &[u8]) -> PrefixEncoder {
        let mut occupancy = [false; 256];
        for &b in text {
            occupancy[b as usize] = true;
        }
        PrefixEncoder::from_occupancy(occupancy)
    }

    #[test]
    fn test_count_fits_the_budget() {
        // 7 bytes: 7^5 = 16807 fits, 7^6 does not
        let dna = encoder(b"$%ACGNT");
        assert_eq!((dna.radix, dna.count, dna.num_buckets()), (7, 5, 16807));

        // 25 bytes: 25^3 = 15625 fits, 25^4 does not
        let protein = encoder(b"$%*-ABCDEFGHIKLMNPQRSTVWY");
        assert_eq!(
            (protein.radix, protein.count, protein.num_buckets()),
            (25, 3, 15625)
        );

        // 256 bytes: 256^2 is the budget exactly
        let bytes: Vec<u8> = (0..=255).collect();
        let bytes = encoder(&bytes);
        assert_eq!(
            (bytes.radix, bytes.count, bytes.num_buckets()),
            (256, 2, 65536)
        );

        // one byte: every key is 0, and count stops at the cap
        let one = encoder(b"$");
        assert_eq!((one.radix, one.count, one.num_buckets()), (1, 16, 1));
    }

    // the keys are written digit by digit, 1 × 6 included
    #[allow(clippy::identity_op)]
    #[test]
    fn test_encode() -> Result<()> {
        // ranks in byte order: $=0 A=1 C=2 G=3 N=4 T=5, so radix 6
        // and 6 characters
        let text = b"ACGTNA$";
        let encoder = encoder(text);
        assert_eq!(encoder.num_buckets(), 46656);

        // A C G T N A = 1 2 3 5 4 1
        let key = ((((1 * 6 + 2) * 6 + 3) * 6 + 5) * 6 + 4) * 6 + 1;
        assert_eq!(encoder.encode(text, 0, None), key);

        // N A $ and three characters past the end: 4 1 0 0 0 0
        let key = (4 * 6 + 1) * 6 * 6 * 6 * 6;
        assert_eq!(encoder.encode(text, 4, None), key);

        // care positions 0, 1 and 3, then three past the mask
        let mask = SeedMask::new("1101")?;
        // A C T = 1 2 5
        let key = ((1 * 6 + 2) * 6 + 5) * 6 * 6 * 6;
        assert_eq!(encoder.encode(text, 0, Some(&mask)), key);
        // A $ and position 8 past the end: 1 0 0
        let key = 6 * 6 * 6 * 6 * 6;
        assert_eq!(encoder.encode(text, 5, Some(&mask)), key);
        Ok(())
    }

    #[test]
    fn test_roll_reaches_the_same_bucket_as_encode() {
        let text = b"ACGTNACGTTTGANNA$";
        let encoder = encoder(text);
        let mut key = encoder.encode(text, 0, None) as usize;
        for i in 0..text.len() {
            assert_eq!(key, encoder.encode(text, i, None) as usize, "at {i}");
            key = encoder.roll(text, i, key);
        }
    }

    #[test]
    fn test_cut_at_cumulative_shares() {
        // 16 suffixes in 3 partitions: a share of 6, so the cuts
        // fall at 6 and 12
        let layout = Layout::cut(&[3, 2, 1, 6, 2, 2], 3);
        assert_eq!(layout.ranges, [0..6, 6..12, 12..16]);
        assert_eq!(layout.buckets, [0..3, 3..4, 4..6]);
        assert_eq!(layout.partition_of, [0, 0, 0, 1, 2, 2]);
        assert_eq!(layout.num_suffixes(), 16);

        // an empty bucket never opens or closes a partition
        let layout = Layout::cut(&[0, 5, 0, 0, 5, 0], 2);
        assert_eq!(layout.ranges, [0..5, 5..10]);
        assert_eq!(layout.buckets, [0..2, 2..5]);
        assert_eq!(layout.partition_of[..5], [0, 0, 1, 1, 1]);

        // more partitions than suffixes
        let layout = Layout::cut(&[1, 0, 1], 1024);
        assert_eq!(layout.ranges, [0..1, 1..2]);
    }

    #[test]
    fn test_counting_sort() {
        let keys = [2, 0, 2, 1, 0];
        let (ordered, starts) = counting_sort(&[10u32, 11, 12, 13, 14], 4, |i| keys[i]);
        assert_eq!(ordered, [11, 14, 13, 10, 12]);
        assert_eq!(starts, [0, 2, 3, 5, 5]);
    }

    #[test]
    fn test_turns() {
        let layout = Layout::cut(&[3, 2, 1, 6, 2, 2], 3);
        assert_eq!(layout.turns(100), vec![0..3; 1]);
        assert_eq!(layout.turns(6), [0..1, 1..2, 2..3]);
        assert_eq!(layout.turns(10), [0..1, 1..3]);
        // a partition over the limit is a turn by itself
        assert_eq!(layout.turns(1), [0..1, 1..2, 2..3]);
    }

    #[test]
    fn test_gather_fills_each_partition() -> Result<()> {
        let text = b"GATTACA$";
        let encoder = encoder(text);
        let keep = [true; 256];

        // 5 ranks and 6 characters: every suffix is its own bucket,
        // in suffix order
        let layout = Layout::new(text, &encoder, &keep, None, 3);
        assert_eq!(layout.ranges, [0..3, 3..6, 6..8]);

        let mut parts = vec![];
        for turn in layout.turns(4) {
            for (mut suffixes, buckets) in
                gather::<u32>(text, &encoder, &keep, None, &layout, turn)?
            {
                assert_eq!(buckets.len(), suffixes.len());
                suffixes.sort();
                parts.push(suffixes);
            }
        }
        // $ A$ ACA$ | ATTACA$ CA$ GATTACA$ | TACA$ TTACA$
        assert_eq!(parts, [vec![4, 6, 7], vec![0, 1, 5], vec![2, 3]]);
        Ok(())
    }
}
