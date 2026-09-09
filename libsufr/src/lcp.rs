//! Functions for LCP computation.

/// Computes the LCP of `a` and `b`, automatically selecting an optimized
/// implementation for the target platform.
pub fn lcp(a: &[u8], b: &[u8]) -> usize {
    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx2") {
            // SAFETY: avx2 is available
            return unsafe { lcp_avx2(a, b) };
        }
    }
    lcp_scalar(a, b)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
/// AVX2 implementation of LCP with an unrolled loop
///
/// SAFETY: This function must only be called after validating that avx2 is actually available
unsafe fn lcp_avx2(a: &[u8], b: &[u8]) -> usize {
    use std::arch::x86_64::*;
    let n = a.len().min(b.len());
    let mut i = 0;
    const UNROLL: usize = 4;
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
    }
    i + lcp_scalar(&a[i..n], &b[i..n])
}

/// Scalar implementation of LCP
pub fn lcp_scalar(a: &[u8], b: &[u8]) -> usize {
    std::iter::zip(a, b).take_while(|(a, b)| a == b).count()
}
