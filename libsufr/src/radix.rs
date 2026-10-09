//! Radix-sort construction of spaced-seed suffix arrays
//!
//! In seed-mask mode the sort key of a suffix is the sequence of text bytes
//! found at the mask's "care" offsets, and suffixes are compared only on
//! those bytes (the LCP is likewise the number of care positions shared).
//! When the key fits in 64 bits it can be packed into a `u64` and the
//! suffix positions sorted by integer key instead of a comparison sort on
//! the text. The same holds for a small maximum query length, where the
//! key is simply the first `max_query_len` bytes.
//!
//! Key layout: every distinct byte in the text gets a rank starting at 1,
//! in byte order; rank 0 is reserved for "past the end of the text" so a
//! truncated suffix sorts before every suffix that extends it. This
//! reproduces the order produced by the merge-sort path, including the
//! tie rule that among equal keys the suffix with the larger text position
//! (the shorter suffix) comes first.
//!
//! Keys are read from a copy of the text converted to ranks in place
//! (`KeyEncoder::convert_to_ranks`). With BMI2 a key is one or two
//! unaligned 8-byte loads and a `pext` each; otherwise a loop over the
//! care offsets.

use anyhow::{bail, Result};
use rayon::prelude::*;
use std::mem;

/// A suffix position together with its packed sort key.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Item<T> {
    pub key: u64,
    pub pos: T,
}

/// Packs the bytes at a set of offsets into a `u64` key.
#[derive(Debug)]
pub(crate) struct KeyEncoder {
    /// Rank of each byte value (0 is reserved for "past end of text")
    rank: [u8; 256],

    /// Bits per symbol in the key
    pub bits: u32,

    /// Total key bits (`bits * positions.len()`)
    pub total_bits: u32,

    /// Offsets of the bytes that make up the key, strictly increasing
    pub positions: Vec<usize>,

    /// One past the largest offset
    pub span: usize,

    /// Length of the text
    pub text_len: usize,

    /// Whether `pext` can be used (BMI2 present and span <= 16)
    use_pext: bool,

    /// `pext` masks for the big-endian 8-byte words at offsets 0 and 8
    pext_mask: [u64; 2],

    /// Key bits produced from the second word
    low_bits: u32,
}

impl KeyEncoder {
    /// Build an encoder for `text` using the given offsets. Fails when the
    /// key would not fit into 64 bits.
    pub fn new(text: &[u8], positions: &[usize]) -> Result<Self> {
        if positions.is_empty() {
            bail!("Radix sort needs at least one key position");
        }

        let present = text
            .par_chunks(1 << 22)
            .map(|chunk| {
                let mut seen = [false; 256];
                for &b in chunk {
                    seen[b as usize] = true;
                }
                seen
            })
            .reduce(
                || [false; 256],
                |mut a, b| {
                    for i in 0..256 {
                        a[i] |= b[i];
                    }
                    a
                },
            );

        let mut rank = [0u8; 256];
        let mut next = 1u32;
        for (b, &seen) in present.iter().enumerate() {
            if seen {
                rank[b] = next as u8;
                next += 1;
            }
        }

        // `next` symbols including the reserved 0: ceil(log2(next)) bits
        let bits = (32 - (next - 1).leading_zeros()).max(1);
        let weight = positions.len();
        let total_bits = bits * weight as u32;
        if total_bits > 64 {
            bail!(
                "Radix sort key needs {weight} symbols x {bits} bits = \
                 {total_bits} bits, which exceeds 64; use a lighter seed mask \
                 or a shorter max query length"
            );
        }
        let span = positions.last().map_or(0, |&p| p + 1);

        // pext masks: a big-endian load puts the byte at offset i in bits
        // 56-8i..64-8i, so selecting the low `bits` bits of each care byte
        // yields the key with offset 0 most significant.
        let symbol_mask = (1u64 << bits) - 1;
        let mut pext_mask = [0u64; 2];
        let mut low_bits = 0;
        for &offset in positions {
            let word = offset / 8;
            let shift = 56 - 8 * (offset % 8) as u32;
            pext_mask[word.min(1)] |= symbol_mask << shift;
            if word == 1 {
                low_bits += bits;
            }
        }
        let use_pext = span <= 16 && has_bmi2();

        Ok(Self {
            rank,
            bits,
            total_bits,
            positions: positions.to_vec(),
            span,
            text_len: text.len(),
            use_pext,
            pext_mask,
            low_bits,
        })
    }

    /// Whether `pext` is used for keys
    pub fn uses_pext(&self) -> bool {
        self.use_pext
    }

    /// Replace every byte of `text` with its rank, in place and in parallel.
    pub fn convert_to_ranks(&self, text: &mut [u8]) {
        text.par_chunks_mut(1 << 22).for_each(|chunk| {
            for b in chunk.iter_mut() {
                *b = self.rank[*b as usize];
            }
        });
    }

    /// Map a predicate on original bytes to a table indexed by rank.
    pub fn rank_table(&self, pred: impl Fn(u8) -> bool) -> [bool; 256] {
        let mut table = [false; 256];
        for b in 0..=255u8 {
            let r = self.rank[b as usize];
            if r != 0 && pred(b) {
                table[r as usize] = true;
            }
        }
        table
    }

    /// The packed key of the suffix starting at `pos`, read from a text
    /// already converted to ranks.
    #[inline(always)]
    pub fn key(&self, ranks: &[u8], pos: usize) -> u64 {
        if self.use_pext && pos + 16 <= self.text_len {
            // Safety: pos + 16 <= text_len == ranks.len()
            unsafe { self.key_pext(ranks, pos) }
        } else if pos + self.span <= self.text_len {
            let mut key = 0u64;
            for &offset in &self.positions {
                // Safety: pos + offset < pos + span <= text_len
                let r = unsafe { *ranks.get_unchecked(pos + offset) };
                key = (key << self.bits) | r as u64;
            }
            key
        } else {
            let mut key = 0u64;
            for &offset in &self.positions {
                let r = ranks.get(pos + offset).copied().unwrap_or(0);
                key = (key << self.bits) | r as u64;
            }
            key
        }
    }

    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "bmi2")]
    #[inline]
    unsafe fn key_pext(&self, ranks: &[u8], pos: usize) -> u64 {
        use std::arch::x86_64::_pext_u64;
        let p = ranks.as_ptr().add(pos);
        let w0 = u64::from_be_bytes(*(p as *const [u8; 8]));
        let hi = _pext_u64(w0, self.pext_mask[0]);
        if self.low_bits == 0 {
            hi
        } else {
            let w1 = u64::from_be_bytes(*(p.add(8) as *const [u8; 8]));
            (hi << self.low_bits) | _pext_u64(w1, self.pext_mask[1])
        }
    }

    #[cfg(not(target_arch = "x86_64"))]
    unsafe fn key_pext(&self, _ranks: &[u8], _pos: usize) -> u64 {
        unreachable!("pext is only enabled on x86_64")
    }

    /// Number of key positions of the suffix at `pos` that fall inside the text.
    #[inline(always)]
    pub fn available(&self, pos: usize) -> usize {
        if pos + self.span <= self.text_len {
            self.positions.len()
        } else {
            self.positions
                .iter()
                .filter(|&&offset| pos + offset < self.text_len)
                .count()
        }
    }

    /// LCP between two suffixes in key symbols (care positions), matching
    /// the definition used by the merge-sort path: the count of leading
    /// equal symbols, never more than either suffix has inside the text.
    #[inline(always)]
    pub fn lcp(&self, a_key: u64, a_pos: usize, b_key: u64, b_pos: usize) -> usize {
        let x = a_key ^ b_key;
        let matched = if x == 0 {
            self.positions.len()
        } else {
            ((x.leading_zeros() - (64 - self.total_bits)) / self.bits) as usize
        };
        matched
            .min(self.available(a_pos))
            .min(self.available(b_pos))
    }
}

#[cfg(target_arch = "x86_64")]
fn has_bmi2() -> bool {
    std::arch::is_x86_feature_detected!("bmi2")
}

#[cfg(not(target_arch = "x86_64"))]
fn has_bmi2() -> bool {
    false
}

/// A raw pointer that may be shared across threads. Callers must write
/// disjoint indices.
struct SharedPtr<T>(*mut T);
unsafe impl<T> Send for SharedPtr<T> {}
unsafe impl<T> Sync for SharedPtr<T> {}
impl<T> SharedPtr<T> {
    #[inline(always)]
    unsafe fn write(&self, index: usize, value: T) {
        self.0.add(index).write(value)
    }
}

/// Number of chunks to split `n` items into for parallel passes.
fn num_chunks(n: usize) -> usize {
    let threads = rayon::current_num_threads();
    (n / (1 << 16)).clamp(1, threads * 4)
}

/// One stable counting-sort pass of `items` into `scratch` by `bucket`,
/// then swap so `items` holds the result. Returns the start index of each
/// bucket (length `num_buckets + 1`). Returns `None` without moving
/// anything when every item falls in the same bucket.
fn counting_pass<T, F>(
    items: &mut Vec<Item<T>>,
    scratch: &mut Vec<Item<T>>,
    num_buckets: usize,
    bucket: F,
) -> Option<Vec<usize>>
where
    T: Copy + Send + Sync,
    F: Fn(&Item<T>) -> usize + Sync,
{
    let n = items.len();
    // Keep the per-chunk histograms small enough to stay in cache: with
    // many buckets, use fewer chunks. Chunks never exceed u32 counts.
    let chunks = num_chunks(n).min((n / (num_buckets * 8)).max(1));
    let chunk_len = n.div_ceil(chunks).max(1);

    // Per-chunk histograms
    let counts: Vec<Vec<u32>> = items
        .par_chunks(chunk_len)
        .map(|chunk| {
            let mut hist = vec![0u32; num_buckets];
            for item in chunk {
                hist[bucket(item)] += 1;
            }
            hist
        })
        .collect();

    // Bucket totals
    let mut totals = vec![0usize; num_buckets];
    for hist in &counts {
        for (t, &c) in totals.iter_mut().zip(hist) {
            *t += c as usize;
        }
    }
    if totals.iter().filter(|&&c| c > 0).count() <= 1 {
        return None;
    }

    // Bucket starts, then per-chunk write offsets (chunk-major within each
    // bucket so the scatter is stable)
    let mut starts = Vec::with_capacity(num_buckets + 1);
    let mut base = 0usize;
    for &t in &totals {
        starts.push(base);
        base += t;
    }
    starts.push(base);
    let mut next = starts[..num_buckets].to_vec();
    let counts: Vec<Vec<usize>> = counts
        .into_iter()
        .map(|hist| {
            hist.into_iter()
                .enumerate()
                .map(|(b, c)| {
                    let offset = next[b];
                    next[b] += c as usize;
                    offset
                })
                .collect()
        })
        .collect();

    // Scatter
    let dst = SharedPtr(scratch.as_mut_ptr());
    items
        .par_chunks(chunk_len)
        .zip(counts.into_par_iter())
        .for_each(|(chunk, mut offsets)| {
            for item in chunk {
                let b = bucket(item);
                // Safety: offsets partition 0..n disjointly across
                // (chunk, bucket) pairs by construction above.
                unsafe { dst.write(offsets[b], *item) };
                offsets[b] += 1;
            }
        });

    mem::swap(items, scratch);
    Some(starts)
}

/// Sort items by key ascending, ties by position descending. Uses the LSD
/// radix sort unless the environment variable `SUFR_RADIX_MSD` is set, in
/// which case the MSD sub-bucket sort is used (kept for comparison: it
/// measured about 30% slower per partition on a human genome partition).
/// LSD is stable, so `items` must already be in descending position order.
pub(crate) fn sort_items<T>(
    items: &mut Vec<Item<T>>,
    scratch: &mut Vec<Item<T>>,
    total_bits: u32,
) where
    T: Copy + Send + Sync + Ord,
{
    if std::env::var_os("SUFR_RADIX_MSD").is_some() {
        msd_sort(items);
    } else {
        lsd_radix_sort(items, scratch, total_bits);
    }
}

/// Stable LSD radix sort of `items` by `key`, considering only the low
/// `total_bits` bits. Uses up to 16-bit digits; a pass whose digit is the
/// same for every item is skipped. `scratch` is resized as needed and may
/// be reused across calls to avoid re-faulting its pages.
pub(crate) fn lsd_radix_sort<T>(
    items: &mut Vec<Item<T>>,
    scratch: &mut Vec<Item<T>>,
    total_bits: u32,
) where
    T: Copy + Send + Sync,
{
    let n = items.len();
    if n < 2 || total_bits == 0 {
        return;
    }

    let max_digit_bits = 16;
    let passes = total_bits.div_ceil(max_digit_bits);
    let digit_bits = total_bits.div_ceil(passes);
    let num_buckets = 1usize << digit_bits;
    let digit_mask = (num_buckets - 1) as u64;
    scratch.resize(n, items[0]);

    for pass in 0..passes {
        let shift = pass * digit_bits;
        counting_pass(items, scratch, num_buckets, |item| {
            ((item.key >> shift) & digit_mask) as usize
        });
    }
}

/// Sort `items` by key ascending, ties by position descending, with one
/// counting pass on the high key bits (up to 2^18 sub-buckets spanning the
/// key range present) followed by an independent sort of each sub-bucket.
pub(crate) fn msd_sort<T>(items: &mut Vec<Item<T>>)
where
    T: Copy + Send + Sync + Ord,
{
    let n = items.len();
    if n < 2 {
        return;
    }
    let by_key_then_pos =
        |a: &Item<T>, b: &Item<T>| a.key.cmp(&b.key).then_with(|| b.pos.cmp(&a.pos));

    // Small inputs: sort directly
    if n <= 1 << 16 {
        items.par_sort_unstable_by(by_key_then_pos);
        return;
    }

    let (min_key, max_key) = items
        .par_iter()
        .map(|item| (item.key, item.key))
        .reduce(|| (u64::MAX, 0), |a, b| (a.0.min(b.0), a.1.max(b.1)));
    if min_key == max_key {
        items.par_sort_unstable_by(by_key_then_pos);
        return;
    }

    // Smallest shift that keeps the number of sub-buckets at or under 2^18
    // (a few hundred to a few thousand items each for a partition of a
    // human genome; the per-chunk histograms then fit in L2)
    let max_buckets = 1u64 << 18;
    let mut shift = 0u32;
    while (max_key >> shift) - (min_key >> shift) >= max_buckets {
        shift += 1;
    }
    let base = min_key >> shift;
    let num_buckets = ((max_key >> shift) - base + 1) as usize;

    let mut scratch: Vec<Item<T>> = vec![items[0]; n];
    let starts = match counting_pass(items, &mut scratch, num_buckets, |item| {
        ((item.key >> shift) - base) as usize
    }) {
        Some(starts) => starts,
        None => {
            items.par_sort_unstable_by(by_key_then_pos);
            return;
        }
    };
    drop(scratch);

    // Sort each sub-bucket independently
    let mut slices: Vec<&mut [Item<T>]> = Vec::with_capacity(num_buckets);
    let mut rest: &mut [Item<T>] = items.as_mut_slice();
    for b in 0..num_buckets {
        let len = starts[b + 1] - starts[b];
        let (head, tail) = rest.split_at_mut(len);
        if len > 1 {
            slices.push(head);
        }
        rest = tail;
    }
    slices.into_par_iter().for_each(|slice| {
        slice.sort_unstable_by(by_key_then_pos);
    });
}

#[cfg(test)]
mod tests {
    use super::{lsd_radix_sort, msd_sort, Item, KeyEncoder};
    use anyhow::Result;
    use pretty_assertions::assert_eq;
    use crate::subsample_tests::Rng;

    fn ranks_of(enc: &KeyEncoder, text: &[u8]) -> Vec<u8> {
        let mut r = text.to_vec();
        enc.convert_to_ranks(&mut r);
        r
    }

    #[test]
    fn test_key_encoder_ranks() -> Result<()> {
        // Distinct bytes: $ A C T -> ranks 1..=4, 3 bits each
        let text = b"ACTA$";
        let enc = KeyEncoder::new(text, &[0, 2])?;
        let ranks = ranks_of(&enc, text);
        assert_eq!(enc.bits, 3);
        assert_eq!(enc.total_bits, 6);
        assert_eq!(enc.span, 3);
        // A=2, T=4
        assert_eq!(enc.key(&ranks, 0), (2 << 3) | 4);
        // A at 3, offset 2 -> 5 is past the end -> 0
        assert_eq!(enc.key(&ranks, 3), 2 << 3);
        assert_eq!(enc.available(3), 1);
        assert_eq!(enc.available(0), 2);
        Ok(())
    }

    #[test]
    fn test_key_encoder_lcp() -> Result<()> {
        let text = b"AAACAAA$";
        let enc = KeyEncoder::new(text, &[0, 2])?;
        let ranks = ranks_of(&enc, text);
        let k = |p| enc.key(&ranks, p);
        // AAA vs ACA: first care position matches, second differs
        assert_eq!(enc.lcp(k(0), 0, k(1), 1), 1);
        // 4: A.A$, 2: A.A -> both symbols match
        assert_eq!(enc.lcp(k(4), 4, k(2), 2), 2);
        // 6: A then past end, 5: A then $ -> only the first can match
        assert_eq!(enc.lcp(k(6), 6, k(5), 5), 1);
        // Two suffixes truncated at the same key position must not count
        // the padding as a match.
        let enc = KeyEncoder::new(text, &[0, 3])?;
        let ranks = ranks_of(&enc, text);
        let k = |p| enc.key(&ranks, p);
        assert_eq!(enc.lcp(k(5), 5, k(6), 6), 1);
        Ok(())
    }

    /// The pext path and the loop must agree on every position for DNA
    /// (3-bit) and protein-like (5-bit) texts and masks spanning up to 16.
    #[test]
    fn test_key_pext_matches_loop() -> Result<()> {
        let mut rng = Rng(21);
        let dna: Vec<u8> = (0..5000)
            .map(|i| if i % 977 == 0 { b'%' } else { b"ACGTN"[rng.below(5)] })
            .chain(std::iter::once(b'$'))
            .collect();
        let prot: Vec<u8> = (0..5000)
            .map(|_| b"ACDEFGHIKLMNPQRSTVWYX%"[rng.below(22)])
            .chain(std::iter::once(b'$'))
            .collect();
        for text in [dna, prot] {
            for mask in ["101", "11011", "11101101101111", "1000000000000001", "111111111111"] {
                let positions: Vec<usize> = mask
                    .bytes()
                    .enumerate()
                    .filter(|(_, b)| *b == b'1')
                    .map(|(i, _)| i)
                    .collect();
                let enc = KeyEncoder::new(&text, &positions)?;
                let ranks = ranks_of(&enc, &text);
                for pos in 0..text.len() {
                    let fast = enc.key(&ranks, pos);
                    let mut slow = 0u64;
                    for &off in &positions {
                        let r = ranks.get(pos + off).copied().unwrap_or(0) as u64;
                        slow = (slow << enc.bits) | r;
                    }
                    assert_eq!(fast, slow, "mask {mask} pos {pos} pext={}", enc.uses_pext());
                }
            }
        }
        Ok(())
    }

    #[test]
    fn test_key_too_wide() {
        let text = b"ACGT$";
        let positions: Vec<usize> = (0..22).collect();
        assert!(KeyEncoder::new(text, &positions).is_err());
    }

    fn random_items(n: u64, key_bits: u32, seed: u64) -> Vec<Item<u64>> {
        let mut rng = Rng(seed);
        let mask = if key_bits == 64 { u64::MAX } else { (1u64 << key_bits) - 1 };
        (0..n)
            .map(|i| Item {
                // Skewed keys so that many ties and many sub-buckets occur
                key: (rng.next() & mask) >> rng.below(3),
                pos: i,
            })
            .collect()
    }

    #[test]
    fn test_lsd_radix_sort_is_stable() {
        let mut items: Vec<Item<u32>> = (0..5000u32)
            .map(|i| Item {
                key: ((i * 7919) % 97) as u64,
                pos: i,
            })
            .collect();
        let mut expected = items.clone();
        expected.sort_by(|a, b| a.key.cmp(&b.key).then(a.pos.cmp(&b.pos)));
        let mut scratch = Vec::new();
        lsd_radix_sort(&mut items, &mut scratch, 7);
        assert_eq!(items, expected);

        // Reuse of a scratch buffer across calls of different sizes
        let mut items = random_items(20000, 44, 1);
        let mut expected = items.clone();
        expected.sort_by(|a, b| a.key.cmp(&b.key).then(a.pos.cmp(&b.pos)));
        let mut scratch = vec![items[0]; 50000];
        lsd_radix_sort(&mut items, &mut scratch, 44);
        assert_eq!(items, expected);
        let mut items = random_items(30000, 44, 2);
        let mut expected = items.clone();
        expected.sort_by(|a, b| a.key.cmp(&b.key).then(a.pos.cmp(&b.pos)));
        lsd_radix_sort(&mut items, &mut scratch, 44);
        assert_eq!(items, expected);
    }

    #[test]
    fn test_msd_sort_orders_by_key_then_descending_pos() {
        for (n, bits) in [(10u64, 12u32), (70000, 12), (300000, 33), (300000, 55), (200000, 64)] {
            let mut items = random_items(n, bits, n);
            let mut expected = items.clone();
            expected.sort_by(|a, b| a.key.cmp(&b.key).then(b.pos.cmp(&a.pos)));
            msd_sort(&mut items);
            assert_eq!(items, expected, "n={n} bits={bits}");
        }
        // All keys equal
        let mut items: Vec<Item<u32>> = (0..100000u32).map(|i| Item { key: 5, pos: i }).collect();
        msd_sort(&mut items);
        assert!(items.windows(2).all(|w| w[0].pos > w[1].pos));
    }
}
