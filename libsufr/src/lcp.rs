//! Functions for LCP computation and caching.

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

/// Computes the LCP of `a` and `b`, automatically selecting an optimized
/// implementation for the target platform.
///
/// If `cache` is provided, it should be provided as `Some((cache, ia, ib))`,
/// `ia` and `ib` are the offsets at which `a` and `b` begin in the original sequence.
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

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
/// AVX2 implementation of LCP with an unrolled loop
///
/// SAFETY: This function must only be called after validating that avx2 is actually available
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
            // SAFETY: the bounds of a and b (256bits = 32 bytes; 4 iterations) are checked at the top of the while loop
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

/// Scalar implementation of LCP
pub fn lcp_scalar(a: &[u8], b: &[u8]) -> usize {
    std::iter::zip(a, b).take_while(|(a, b)| a == b).count()
}

#[cfg(test)]
mod tests {
    use super::*;

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
}
