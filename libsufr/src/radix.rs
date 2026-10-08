//! Radix-sort construction of spaced-seed suffix arrays
//!
//! In seed-mask mode the sort key of a suffix is the sequence of text bytes
//! found at the mask's "care" offsets, and suffixes are compared only on
//! those bytes (the LCP is likewise the number of care positions shared).
//! When the key fits in 64 bits it can be packed into a `u64` and the
//! suffix positions sorted with an LSD radix sort instead of a comparison
//! sort. The same holds for a small maximum query length, where the key is
//! simply the first `max_query_len` bytes.
//!
//! Key layout: every distinct byte in the text gets a rank starting at 1,
//! in byte order; rank 0 is reserved for "past the end of the text" so a
//! truncated suffix sorts before every suffix that extends it. This
//! reproduces the order produced by the merge-sort path, including the
//! tie rule that among equal keys the suffix with the larger text position
//! (the shorter suffix) comes first.

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

        Ok(Self {
            rank,
            bits,
            total_bits,
            positions: positions.to_vec(),
            span: positions.last().map_or(0, |&p| p + 1),
            text_len: text.len(),
        })
    }

    /// The packed key of the suffix starting at `pos`.
    #[inline(always)]
    pub fn key(&self, text: &[u8], pos: usize) -> u64 {
        let mut key = 0u64;
        if pos + self.span <= self.text_len {
            for &offset in &self.positions {
                // Safety: pos + offset < pos + span <= text_len
                let b = unsafe { *text.get_unchecked(pos + offset) };
                key = (key << self.bits) | self.rank[b as usize] as u64;
            }
        } else {
            for &offset in &self.positions {
                let r = if pos + offset < self.text_len {
                    self.rank[text[pos + offset] as usize] as u64
                } else {
                    0
                };
                key = (key << self.bits) | r;
            }
        }
        key
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

/// Stable LSD radix sort of `items` by `key`, considering only the low
/// `total_bits` bits. Uses up to 16-bit digits; a pass whose digit is the
/// same for every item is skipped.
pub(crate) fn lsd_radix_sort<T>(items: &mut Vec<Item<T>>, total_bits: u32)
where
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

    let chunks = num_chunks(n);
    let chunk_len = n.div_ceil(chunks);

    let mut scratch: Vec<Item<T>> = vec![items[0]; n];

    for pass in 0..passes {
        let shift = pass * digit_bits;
        let digit = |item: &Item<T>| ((item.key >> shift) & digit_mask) as usize;

        // Per-chunk histograms
        let mut counts: Vec<Vec<usize>> = items
            .par_chunks(chunk_len)
            .map(|chunk| {
                let mut hist = vec![0usize; num_buckets];
                for item in chunk {
                    hist[digit(item)] += 1;
                }
                hist
            })
            .collect();

        // Bucket totals; skip the pass when every item has the same digit
        let mut totals = vec![0usize; num_buckets];
        for hist in &counts {
            for (t, &c) in totals.iter_mut().zip(hist) {
                *t += c;
            }
        }
        if totals.iter().filter(|&&c| c > 0).count() <= 1 {
            continue;
        }

        // Turn counts into starting offsets, chunk-major within each bucket
        // so the scatter is stable.
        let mut base = 0usize;
        let mut starts = vec![0usize; num_buckets];
        for (b, &t) in totals.iter().enumerate() {
            starts[b] = base;
            base += t;
        }
        for hist in counts.iter_mut() {
            for (b, c) in hist.iter_mut().enumerate() {
                let count = *c;
                *c = starts[b];
                starts[b] += count;
            }
        }

        // Scatter
        let dst = SharedPtr(scratch.as_mut_ptr());
        items
            .par_chunks(chunk_len)
            .zip(counts.into_par_iter())
            .for_each(|(chunk, mut offsets)| {
                for item in chunk {
                    let b = digit(item);
                    // Safety: offsets partition 0..n disjointly across
                    // (chunk, bucket) pairs by construction above.
                    unsafe { dst.write(offsets[b], *item) };
                    offsets[b] += 1;
                }
            });

        mem::swap(items, &mut scratch);
    }
}

#[cfg(test)]
mod tests {
    use super::{lsd_radix_sort, Item, KeyEncoder};
    use anyhow::Result;
    use pretty_assertions::assert_eq;

    #[test]
    fn test_key_encoder_ranks() -> Result<()> {
        // Distinct bytes: $ A C T -> ranks 1..=4, 3 bits each
        let text = b"ACTA$";
        let enc = KeyEncoder::new(text, &[0, 2])?;
        assert_eq!(enc.bits, 3);
        assert_eq!(enc.total_bits, 6);
        assert_eq!(enc.span, 3);
        // A=2, T=4
        assert_eq!(enc.key(text, 0), (2 << 3) | 4);
        // A at 3, offset 2 -> 5 is past the end -> 0
        assert_eq!(enc.key(text, 3), 2 << 3);
        assert_eq!(enc.available(3), 1);
        assert_eq!(enc.available(0), 2);
        Ok(())
    }

    #[test]
    fn test_key_encoder_lcp() -> Result<()> {
        let text = b"AAACAAA$";
        let enc = KeyEncoder::new(text, &[0, 2])?;
        let k = |p| enc.key(text, p);
        // AAA vs ACA: first care position matches, second differs
        assert_eq!(enc.lcp(k(0), 0, k(1), 1), 1);
        // 4: A.A$, 2: A.A -> both symbols match
        assert_eq!(enc.lcp(k(4), 4, k(2), 2), 2);
        // 6: A then past end, 5: A then $ -> only the first can match
        assert_eq!(enc.lcp(k(6), 6, k(5), 5), 1);
        // Two suffixes truncated at the same key position must not count
        // the padding as a match.
        let enc = KeyEncoder::new(text, &[0, 3])?;
        let k = |p| enc.key(text, p);
        assert_eq!(enc.lcp(k(5), 5, k(6), 6), 1);
        Ok(())
    }

    #[test]
    fn test_key_too_wide() {
        let text = b"ACGT$";
        let positions: Vec<usize> = (0..22).collect();
        assert!(KeyEncoder::new(text, &positions).is_err());
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
        lsd_radix_sort(&mut items, 7);
        assert_eq!(items, expected);

        // Wider keys, several passes
        let mut items: Vec<Item<u64>> = (0..20000u64)
            .map(|i| Item {
                key: (i.wrapping_mul(0x9E3779B97F4A7C15)) >> 20,
                pos: i,
            })
            .collect();
        let mut expected = items.clone();
        expected.sort_by(|a, b| a.key.cmp(&b.key).then(a.pos.cmp(&b.pos)));
        lsd_radix_sort(&mut items, 44);
        assert_eq!(items, expected);
    }
}
