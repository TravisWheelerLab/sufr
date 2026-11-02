//! Functions for LCP computation and caching.

use crate::format::SeedMask;
use std::cmp::min;

/// A cache of LCP values for optimized calculation.
///
/// Given two indices i and j, with an LCP of l at that position,
/// the LCP of i+1 and j+1 is l-1, of i+2 and j+2 is l-2, and so on.
/// This structure caches LCP values using the difference between i and j
/// as a key, so that LCP values can be reused.
pub struct LcpCache {
    // (2^CACHE_BITS) entries
    // (distance between two starts, smaller start value, lcp)
    cache: Box<[(usize, usize, usize)]>,
}

impl Default for LcpCache {
    fn default() -> Self {
        LcpCache {
            cache: vec![(usize::MAX, 0, 0); Self::CACHE_SIZE].into_boxed_slice(),
        }
    }
}

impl LcpCache {
    /// Index at which to probe an LCP cache, if given.
    /// Must be a multiple of the (unrolled) loop size below.
    pub const PROBE_AT: usize = 1024;
    /// Minimum LCP that should be cached. Must be >= PROBE_AT.
    pub const MIN_CACHEABLE: usize = 1024;

    const CACHE_BITS: usize = 16;
    const CACHE_SIZE: usize = 1 << Self::CACHE_BITS;

    fn idx(d: usize) -> usize {
        // Fibonacci hash for 64-bit integers
        ((d as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15) >> (64 - Self::CACHE_BITS))
            as usize
    }

    /// Adds the LCP of i,j to the cache for later lookup.
    /// Note that this must be the *exact* LCP between i and j,
    /// i.e. not capped to a max length and including any prior skip length.
    pub fn set(&mut self, i: usize, j: usize, lcp: usize) {
        let first = i.min(j);
        let d = i.abs_diff(j);
        let idx = Self::idx(d);
        let loc = &mut self.cache[idx];

        // Eviction conditions:
        // * longer length
        // * both i and j are after the currently-cached area
        // * different distance between
        if lcp >= loc.2 || d != loc.0 || first > loc.1 + loc.2 {
            *loc = (d, first, lcp);
        }
    }

    /// Attempts to find the LCP of i,j using the cache.
    /// Returns `Some(lcp)` if it can be computed from the cache, otherwise `None`.
    pub fn get(&self, i: usize, j: usize) -> Option<usize> {
        let first = i.min(j);
        let d = i.abs_diff(j);
        let idx = Self::idx(d);
        let (loc_d, loc_start, loc_lcp) = self.cache[idx];
        if d == loc_d && first >= loc_start && first <= loc_start + loc_lcp {
            Some(loc_lcp - (first - loc_start))
        } else {
            None
        }
    }
}

const _: () = assert!(LcpCache::MIN_CACHEABLE >= LcpCache::PROBE_AT);

/// The LCP of the suffixes of `text` at `start1` and `start2`, which share
/// their first `skip` bytes, comparing at most `len` bytes past `skip`.
#[inline(always)]
pub fn suffix_lcp(
    text: &[u8],
    start1: usize,
    start2: usize,
    skip: usize,
    len: usize,
    cache: Option<&mut LcpCache>,
) -> usize {
    let text_len = text.len();
    let from1 = start1 + skip;
    let from2 = start2 + skip;
    let end1 = min(from1 + len, text_len);
    let end2 = min(from2 + len, text_len);
    let lcp = lcp(
        &text[from1..end1],
        &text[from2..end2],
        cache.as_deref().map(|c| (c, from1, from2)),
    );

    // a capped LCP is not the true LCP of the two
    // suffixes, and the cache holds only true LCPs
    if skip + lcp >= LcpCache::MIN_CACHEABLE && lcp < len {
        if let Some(c) = cache {
            c.set(start1, start2, skip + lcp);
        }
    }
    skip + lcp
}

/// The number of leading `positions` at which `a` and `b` hold the same
/// byte, stopping at the first position past the end of either.
#[inline(always)]
pub fn masked_lcp(a: &[u8], b: &[u8], positions: &[usize]) -> usize {
    positions
        .iter()
        .take_while(|&&offset| {
            matches!((a.get(offset), b.get(offset)), (Some(x), Some(y)) if x == y)
        })
        .count()
}

/// The offset of the byte that follows a masked LCP of `lcp` care
/// positions: the next care position, or one past the last one.
pub fn mask_offset(lcp: usize, mask: &SeedMask) -> usize {
    if lcp == 0 || lcp > mask.bytes.len() {
        return lcp;
    }
    match mask.positions.get(lcp) {
        Some(&next) => next,
        None => mask.positions[lcp - 1] + 1,
    }
}

/// The LCP of `a` and `b`, using AVX2 when the CPU has it. `cache` holds
/// the cache and the text positions where `a` and `b` start.
pub fn lcp(a: &[u8], b: &[u8], cache: Option<(&LcpCache, usize, usize)>) -> usize {
    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx2") {
            // SAFETY: avx2 is available
            return unsafe { lcp_avx2(a, b, cache) };
        }
    }
    // TODO: aarch64 implementation (including LCP cache)

    lcp_scalar(a, b)
}

/// The LCP of `a` and `b`, comparing 32 bytes at a time.
///
/// # Safety
///
/// The CPU must support AVX2.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn lcp_avx2(
    a: &[u8],
    b: &[u8],
    cache: Option<(&LcpCache, usize, usize)>,
) -> usize {
    use std::arch::x86_64::*;
    let n = a.len().min(b.len());
    let mut i = 0;

    const UNROLL: usize = 4;
    const _: () = assert!(LcpCache::PROBE_AT % (UNROLL * 32) == 0);

    while i + UNROLL * 32 <= n {
        for j in 0..UNROLL {
            let off = i + j * 32;
            // SAFETY: the loop condition keeps the 32 bytes at
            //         `off` within both `a` and `b`
            let va = unsafe { _mm256_loadu_si256(a.as_ptr().add(off) as *const _) };
            let vb = unsafe { _mm256_loadu_si256(b.as_ptr().add(off) as *const _) };
            let eq = _mm256_cmpeq_epi8(va, vb);
            let mask = _mm256_movemask_epi8(eq) as u32;
            if mask != 0xFFFF_FFFF {
                return off + (!mask).trailing_zeros() as usize;
            }
        }
        i += UNROLL * 32;
        if i % LcpCache::PROBE_AT == 0 {
            if let Some((lcp_cache, ia, ib)) = cache {
                if let Some(lcp) = lcp_cache.get(ia + i, ib + i) {
                    return (i + lcp).min(n);
                }
            }
        }
    }
    i + lcp_scalar(&a[i..n], &b[i..n])
}

/// The LCP of `a` and `b`, comparing a byte at a time.
pub fn lcp_scalar(a: &[u8], b: &[u8]) -> usize {
    std::iter::zip(a, b).take_while(|(a, b)| a == b).count()
}

#[cfg(test)]
mod tests {
    use super::*;
    use anyhow::Result;

    #[test]
    fn test_lcp_cache() {
        let mut cache = LcpCache::default();
        cache.set(4, 2000, 1000);
        cache.set(8, 4000, 1000);
        cache.set(12, 8000, 1000);

        assert_eq!(cache.get(10, 2006), Some(994));
        assert_eq!(cache.get(10, 4002), Some(998));
        assert_eq!(cache.get(20, 8008), Some(992));
        assert_eq!(cache.get(20, 3000), None);
    }

    #[test]
    #[cfg(target_arch = "x86_64")] // TODO: consider implementing the cache on other architectures
    fn test_lcp_with_cache() {
        let mut cache = LcpCache::default();
        let mut text = vec![b'A'; 4096];
        text[2047] = b'C';

        // Set a "lie" in the cache to confirm that it is used
        cache.set(0, 2048, 1024);
        assert_eq!(2047, lcp(&text[0..], &text[2048..], None));
        assert_eq!(
            1024,
            lcp(&text[0..], &text[2048..], Some((&cache, 0, 2048)))
        );
    }

    #[test]
    fn test_suffix_lcp() {
        let text = b"ACGTNNACGT$";
        assert_eq!(suffix_lcp(text, 0, 6, 0, text.len(), None), 4);
        assert_eq!(suffix_lcp(text, 0, 6, 0, 2, None), 2);
        assert_eq!(suffix_lcp(text, 0, 6, 2, 1, None), 3);
    }

    #[test]
    fn test_mask_offset() -> Result<()> {
        let seed_mask = SeedMask::new("101")?;
        assert_eq!(mask_offset(0, &seed_mask), 0);
        assert_eq!(mask_offset(1, &seed_mask), 2);
        assert_eq!(mask_offset(2, &seed_mask), 3);

        let seed_mask = SeedMask::new("11011")?;
        assert_eq!(mask_offset(0, &seed_mask), 0);
        assert_eq!(mask_offset(1, &seed_mask), 1);
        assert_eq!(mask_offset(2, &seed_mask), 3);
        assert_eq!(mask_offset(3, &seed_mask), 4);
        assert_eq!(mask_offset(4, &seed_mask), 5);

        let seed_mask = SeedMask::new("10011001")?;
        assert_eq!(mask_offset(0, &seed_mask), 0);
        assert_eq!(mask_offset(1, &seed_mask), 3);
        assert_eq!(mask_offset(2, &seed_mask), 4);
        assert_eq!(mask_offset(3, &seed_mask), 7);
        assert_eq!(mask_offset(4, &seed_mask), 8);

        Ok(())
    }
}
