use crate::{format::SeedMask, io::Int};
use anyhow::{anyhow, Result};
use log::info;
use rayon::prelude::*;
use std::{ops::Range, sync::Mutex, time::Instant};

/// The suffixes one task of a pass over the text takes.
const CHUNK_SIZE: usize = 1 << 19;

/// The rank of each byte of the text in byte order, packed into a key:
/// the ranks of a suffix's first `key_chars` characters at `bits` each
/// from the top of a `u64`, so keys compare as the characters do and
/// two keys share `leading_zeros(a ^ b) / bits` characters.
pub(super) struct PrefixEncoder {
    table: [u8; 256],

    /// The bits per rank.
    bits: u32,

    /// The characters in a key.
    key_chars: usize,

    /// The characters in a bucket number, the key's top ones.
    bucket_chars: usize,

    /// The unused low bits of a key.
    shift: u32,

    /// Keeps the first `max_query_len` characters of a key, or all.
    cap: u64,
}

impl PrefixEncoder {
    /// The encoder for a text whose bytes `occupancy` marks, sorted by
    /// its first `max_query_len` characters, or all of them for 0.
    pub(super) fn from_occupancy(occupancy: [bool; 256], max_query_len: usize) -> Self {
        let mut table = [0u8; 256];
        let mut radix: usize = 0;
        for (c, &seen) in occupancy.iter().enumerate() {
            if seen {
                table[c] = radix as u8;
                radix += 1;
            }
        }
        let bits = radix.next_power_of_two().trailing_zeros().max(1);
        let key_chars = 64 / bits as usize;
        let bucket_chars = 16 / bits as usize;
        let shift = 64 - bits * key_chars as u32;
        let cap = if max_query_len > 0 && max_query_len < key_chars {
            u64::MAX << (64 - bits * max_query_len as u32)
        } else {
            u64::MAX
        };
        info!(
            "{radix} distinct bytes at {bits} bits: {key_chars}-character keys, \
             {bucket_chars}-character buckets"
        );
        Self {
            table,
            bits,
            key_chars,
            bucket_chars,
            shift,
            cap,
        }
    }

    /// The number of distinct bucket numbers.
    pub(super) fn num_buckets(&self) -> usize {
        1 << (self.bits * self.bucket_chars as u32)
    }

    /// `key` with the characters past `max_query_len` cleared.
    #[inline(always)]
    pub(super) fn capped(&self, key: u64) -> u64 {
        key & self.cap
    }

    /// The bucket number of a suffix with the capped key `key`.
    #[inline(always)]
    pub(super) fn bucket(&self, key: u64) -> usize {
        (key >> (64 - self.bits * self.bucket_chars as u32)) as usize
    }

    /// The characters two suffixes with the unequal keys `a` and `b`
    /// share.
    #[inline(always)]
    pub(super) fn shared(&self, a: u64, b: u64) -> usize {
        ((a ^ b).leading_zeros() / self.bits) as usize
    }

    /// The key of the suffix at `i` before the cap: its first
    /// `key_chars` characters, or those at the care positions of `mask`.
    pub(super) fn encode(&self, text: &[u8], i: usize, mask: Option<&SeedMask>) -> u64 {
        let mut key = 0u64;
        for j in 0..self.key_chars {
            // a character past the end of the text or the mask
            // ranks 0, so a short suffix keeps every digit
            let c = match mask {
                Some(mask) => mask.positions.get(j).and_then(|&off| text.get(i + off)),
                None => text.get(i + j),
            };
            let rank = c.map_or(0, |&c| self.table[c as usize] as u64);
            key = (key << self.bits) | rank;
        }
        key << self.shift
    }

    /// The key of the suffix after the one at `i` with the key `key`
    /// before the cap, without a mask: the first character shifts out
    /// of the top and the character `key_chars` ahead comes in at the
    /// bottom.
    #[inline(always)]
    pub(super) fn roll(&self, text: &[u8], i: usize, key: u64) -> u64 {
        let incoming = text
            .get(i + self.key_chars)
            .map_or(0, |&c| self.table[c as usize] as u64);
        (key << self.bits) | (incoming << self.shift)
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
                                    counts[encoder.bucket(encoder.capped(key))] += 1;
                                }
                            }
                        }
                        None => {
                            let mut key = encoder.encode(text, start, None);
                            for i in start..end {
                                if keep[text[i] as usize] {
                                    counts[encoder.bucket(encoder.capped(key))] += 1;
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

/// The suffixes of each partition of `partitions` with their keys, from
/// one scan of `text` keeping those whose first byte `keep` marks.
pub(super) fn gather<W: Int>(
    text: &[u8],
    encoder: &PrefixEncoder,
    keep: &[bool; 256],
    mask: Option<&SeedMask>,
    layout: &Layout,
    partitions: Range<usize>,
) -> Result<Vec<(Vec<W>, Vec<u64>)>> {
    let now = Instant::now();
    let first = partitions.start;
    let count = partitions.len();

    // one compare for "in this turn's buckets", so the branch is
    // rarely taken whichever turn this is; as lo <= b && b < hi
    // the first test is a coin flip in the middle turns and the
    // scan ran three times slower there
    let lo = layout.buckets[first].start;
    let width = layout.buckets[partitions.end - 1].end - lo;
    let gathered: Vec<Mutex<(Vec<W>, Vec<u64>)>> = partitions
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
            || (vec![W::default(); CHUNK_SIZE], vec![0u64; CHUNK_SIZE]),
            |(suffixes, keys), chunk| -> Result<()> {
                let start = chunk * CHUNK_SIZE;
                let end = (start + CHUNK_SIZE).min(text.len());
                let mut n = 0;
                match mask {
                    Some(mask) => {
                        for i in start..end {
                            let key =
                                encoder.capped(encoder.encode(text, i, Some(mask)));
                            let bucket = encoder.bucket(key);
                            if keep[text[i] as usize] && bucket.wrapping_sub(lo) < width
                            {
                                suffixes[n] = W::from_usize(i);
                                keys[n] = key;
                                n += 1;
                            }
                        }
                    }
                    None => {
                        let mut key = encoder.encode(text, start, None);
                        for i in start..end {
                            let capped = encoder.capped(key);
                            let bucket = encoder.bucket(capped);
                            if keep[text[i] as usize] && bucket.wrapping_sub(lo) < width
                            {
                                suffixes[n] = W::from_usize(i);
                                keys[n] = capped;
                                n += 1;
                            }
                            key = encoder.roll(text, i, key);
                        }
                    }
                }
                let suffixes = &suffixes[..n];
                let keys = &keys[..n];

                let partition =
                    |i: usize| layout.partition_of[encoder.bucket(keys[i])] - first;
                let (suffixes, starts) = counting_sort(suffixes, count, partition);
                let (keys, _) = counting_sort(keys, count, partition);
                for (p, run) in starts.windows(2).enumerate() {
                    let run = run[0]..run[1];
                    if !run.is_empty() {
                        let mut g = gathered[p].lock().map_err(|e| anyhow!("{e}"))?;
                        g.0.extend_from_slice(&suffixes[run.clone()]);
                        g.1.extend_from_slice(&keys[run]);
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

    fn encoder(text: &[u8], max_query_len: usize) -> PrefixEncoder {
        let mut occupancy = [false; 256];
        for &b in text {
            occupancy[b as usize] = true;
        }
        PrefixEncoder::from_occupancy(occupancy, max_query_len)
    }

    /// `ranks` packed from the top of a u64 at `bits` each.
    fn packed(ranks: &[u64], bits: u32) -> u64 {
        ranks
            .iter()
            .enumerate()
            .map(|(j, &r)| r << (64 - bits * (j as u32 + 1)))
            .sum()
    }

    #[test]
    fn test_widths() {
        let dna = encoder(b"$%ACGNT", 0);
        assert_eq!((dna.bits, dna.key_chars, dna.bucket_chars), (3, 21, 5));
        assert_eq!(dna.num_buckets(), 32768);

        let protein = encoder(b"$%*-ABCDEFGHIKLMNPQRSTVWY", 0);
        assert_eq!(
            (protein.bits, protein.key_chars, protein.bucket_chars),
            (5, 12, 3)
        );

        let bytes: Vec<u8> = (0..=255).collect();
        let bytes = encoder(&bytes, 0);
        assert_eq!((bytes.bits, bytes.key_chars, bytes.bucket_chars), (8, 8, 2));
        assert_eq!(bytes.num_buckets(), 65536);

        let one = encoder(b"$", 0);
        assert_eq!((one.bits, one.key_chars, one.bucket_chars), (1, 64, 16));
    }

    #[test]
    fn test_encode() -> Result<()> {
        // ranks in byte order: $=0 A=1 C=2 G=3 N=4 T=5, 3 bits each
        let text = b"ACGTNA$";
        let encoder = encoder(text, 0);

        // A C G T N A $ and 14 characters past the end
        let key = encoder.encode(text, 0, None);
        assert_eq!(key, packed(&[1, 2, 3, 5, 4, 1, 0], 3));
        assert_eq!(
            encoder.bucket(key),
            packed(&[1, 2, 3, 5, 4], 3) as usize >> 49
        );

        // N A $ and the rest past the end
        assert_eq!(encoder.encode(text, 4, None), packed(&[4, 1, 0], 3));

        // care positions 0, 1 and 3, then past the mask: A C T
        let mask = SeedMask::new("1101")?;
        assert_eq!(encoder.encode(text, 0, Some(&mask)), packed(&[1, 2, 5], 3));
        // A $ and position 8 past the end
        assert_eq!(encoder.encode(text, 5, Some(&mask)), packed(&[1], 3));

        // a cap keeps that many characters
        let capped = super::test::encoder(text, 2);
        assert_eq!(
            capped.capped(capped.encode(text, 0, None)),
            packed(&[1, 2], 3)
        );
        Ok(())
    }

    #[test]
    fn test_roll_reaches_the_same_key_as_encode() {
        let text = b"ACGTNACGTTTGANNACGTACGTACGTAAA$";
        for max_query_len in [0, 4] {
            let encoder = encoder(text, max_query_len);
            let mut key = encoder.encode(text, 0, None);
            for i in 0..text.len() {
                assert_eq!(key, encoder.encode(text, i, None), "at {i}");
                key = encoder.roll(text, i, key);
            }
        }
    }

    #[test]
    fn test_shared_counts_equal_leading_characters() {
        let text = b"ACGTACGA$";
        let encoder = encoder(text, 0);
        let a = encoder.encode(text, 0, None);
        let b = encoder.encode(text, 4, None);
        // ACGTACGA$ and ACGA$ share ACG
        assert_eq!(encoder.shared(a, b), 3);
        // the lowest bit of a character deciding is still that character
        let c = packed(&[1, 2, 3, 4], 3);
        let d = packed(&[1, 2, 3, 5], 3);
        assert_eq!(encoder.shared(c, d), 3);
        assert_eq!(encoder.shared(a, encoder.encode(text, 8, None)), 0);
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
    fn test_turns() {
        let layout = Layout::cut(&[3, 2, 1, 6, 2, 2], 3);
        assert_eq!(layout.turns(100), vec![0..3; 1]);
        assert_eq!(layout.turns(6), [0..1, 1..2, 2..3]);
        assert_eq!(layout.turns(10), [0..1, 1..3]);
        // a partition over the limit is a turn by itself
        assert_eq!(layout.turns(1), [0..1, 1..2, 2..3]);
    }

    #[test]
    fn test_counting_sort() {
        let keys = [2, 0, 2, 1, 0];
        let (ordered, starts) = counting_sort(&[10u32, 11, 12, 13, 14], 4, |i| keys[i]);
        assert_eq!(ordered, [11, 14, 13, 10, 12]);
        assert_eq!(starts, [0, 2, 3, 5, 5]);
    }

    #[test]
    fn test_gather_fills_each_partition() -> Result<()> {
        let text = b"GATTACA$";
        let encoder = encoder(text, 0);
        let keep = [true; 256];

        // 5 ranks at 3 bits, 5-character buckets: every suffix is its
        // own bucket, in suffix order
        let layout = Layout::new(text, &encoder, &keep, None, 3);
        assert_eq!(layout.ranges, [0..3, 3..6, 6..8]);

        let mut parts = vec![];
        for turn in layout.turns(4) {
            for (mut suffixes, keys) in
                gather::<u32>(text, &encoder, &keep, None, &layout, turn)?
            {
                for (&s, &k) in suffixes.iter().zip(&keys) {
                    assert_eq!(
                        k,
                        encoder.capped(encoder.encode(text, s as usize, None))
                    );
                }
                suffixes.sort();
                parts.push(suffixes);
            }
        }
        // $ A$ ACA$ | ATTACA$ CA$ GATTACA$ | TACA$ TTACA$
        assert_eq!(parts, [vec![4, 6, 7], vec![0, 1, 5], vec![2, 3]]);
        Ok(())
    }
}
