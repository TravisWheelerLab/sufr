use super::{partition::Partition, scratch::Scratch};
use crate::{
    format::SeedMask,
    io::Int,
    lcp::{mask_offset, masked_lcp, suffix_lcp, LcpCache},
};
use anyhow::Result;
use log::info;
use rayon::prelude::*;
use std::{
    cell::RefCell,
    cmp::{max, min, Ordering},
    mem,
    ops::Range,
    time::Instant,
};
use thread_local::ThreadLocal;

/// Sort the suffixes of each group of buffers into one partition, in
/// order, dropping the groups that hold none.
pub(super) fn sort_partitions<T: Int>(
    text: &[u8],
    groups: Vec<Vec<Scratch<T>>>,
    mask: Option<&SeedMask>,
    max_query_len: usize,
    n_ranges: &[Range<usize>],
    write_lcp: bool,
) -> Result<Vec<Partition<T>>> {
    let num_partitions = groups.len();
    let total_sort_time = Instant::now();
    let mut partitions: Vec<Option<Partition<T>>> =
        (0..num_partitions).map(|_| None).collect();

    let lcp_cache: ThreadLocal<RefCell<LcpCache>> = ThreadLocal::new();
    partitions
        .par_iter_mut()
        .zip(groups)
        .enumerate()
        .try_for_each(|(order, (partition, buffers))| -> Result<()> {
            *partition = sort_partition(
                text,
                order,
                buffers,
                mask,
                max_query_len,
                n_ranges,
                write_lcp,
                &lcp_cache,
            )?;
            Ok(())
        })?;

    // drop the empty partitions and put the rest in order
    let mut partitions: Vec<_> = partitions.into_iter().flatten().collect();
    partitions.sort_by_key(|p| p.order);

    let sizes: Vec<_> = partitions.iter().map(|p| p.len).collect();
    let total_size = sizes.iter().sum::<usize>();
    info!(
        "Sorted {total_size} suffixes in {num_partitions} partitions (avg {}) in {:?}",
        total_size / num_partitions,
        total_sort_time.elapsed()
    );
    Ok(partitions)
}

/// Sort the suffixes in `buffers` into the partition at `order`, or
/// `None` when there are none.
#[allow(clippy::too_many_arguments)]
fn sort_partition<T: Int>(
    text: &[u8],
    order: usize,
    buffers: Vec<Scratch<T>>,
    mask: Option<&SeedMask>,
    max_query_len: usize,
    n_ranges: &[Range<usize>],
    write_lcp: bool,
    lcp_cache: &ThreadLocal<RefCell<LcpCache>>,
) -> Result<Option<Partition<T>>> {
    let Some(kind) = buffers.first().map(Scratch::kind) else {
        return Ok(None);
    };
    let mut part_sa = vec![];
    for buffer in buffers {
        part_sa.extend_from_slice(&buffer.consume()?);
    }

    let len = part_sa.len();
    if len == 0 {
        return Ok(None);
    }

    let mut sa_w = part_sa.clone();
    let mut lcp = vec![T::default(); len];
    let mut lcp_w = vec![T::default(); len];
    merge_sort(
        text,
        Arrays {
            sa: &mut sa_w,
            lcp: &mut lcp_w,
        },
        Arrays {
            sa: &mut part_sa,
            lcp: &mut lcp,
        },
        mask,
        max_query_len,
        n_ranges,
        lcp_cache,
    );

    let mut sa_buffer = Scratch::new(kind);
    sa_buffer.extend_from_slice(&part_sa)?;
    let lcp_buffer = if write_lcp {
        let mut lcp_buffer = Scratch::new(kind);
        lcp_buffer.extend_from_slice(&lcp)?;
        Some(lcp_buffer)
    } else {
        None
    };

    Ok(Some(Partition {
        order,
        len,
        first_suffix: part_sa.first().unwrap().to_usize(),
        last_suffix: part_sa.last().unwrap().to_usize(),
        sa_buffer,
        lcp_buffer,
    }))
}

/// A slice of the suffix array and the LCPs beside it.
struct Arrays<'a, T> {
    sa: &'a mut [T],
    lcp: &'a mut [T],
}

impl<T> Arrays<'_, T> {
    fn split_at_mut(&mut self, mid: usize) -> (Arrays<'_, T>, Arrays<'_, T>) {
        let (sa_l, sa_r) = self.sa.split_at_mut(mid);
        let (lcp_l, lcp_r) = self.lcp.split_at_mut(mid);
        (
            Arrays {
                sa: sa_l,
                lcp: lcp_l,
            },
            Arrays {
                sa: sa_r,
                lcp: lcp_r,
            },
        )
    }
}

/// Sort the suffixes in `out` with their LCPs, using `work`, which
/// starts with the same suffixes, as scratch space.
fn merge_sort<T: Int>(
    text: &[u8],
    mut work: Arrays<'_, T>,
    mut out: Arrays<'_, T>,
    mask: Option<&SeedMask>,
    max_query_len: usize,
    n_ranges: &[Range<usize>],
    lcp_cache: &ThreadLocal<RefCell<LcpCache>>,
) {
    let n = work.sa.len();
    if n == 1 {
        out.lcp[0] = T::default();
    } else {
        const NESTED_PAR_GRAIN_SIZE: usize = 1 << 13;
        let mid = n / 2;
        let (work_l, work_r) = work.split_at_mut(mid);
        let (out_l, out_r) = out.split_at_mut(mid);

        // each half is sorted into `work`, swapping the
        // roles, so the merge below reads `work`
        if mid < NESTED_PAR_GRAIN_SIZE {
            merge_sort(
                text,
                out_l,
                work_l,
                mask,
                max_query_len,
                n_ranges,
                lcp_cache,
            );
            merge_sort(
                text,
                out_r,
                work_r,
                mask,
                max_query_len,
                n_ranges,
                lcp_cache,
            );
        } else {
            rayon::join(
                || {
                    merge_sort(
                        text,
                        out_l,
                        work_l,
                        mask,
                        max_query_len,
                        n_ranges,
                        lcp_cache,
                    )
                },
                || {
                    merge_sort(
                        text,
                        out_r,
                        work_r,
                        mask,
                        max_query_len,
                        n_ranges,
                        lcp_cache,
                    )
                },
            );
        }

        merge(
            text,
            work,
            mid,
            out,
            mask,
            max_query_len,
            n_ranges,
            &mut lcp_cache.get_or_default().borrow_mut(),
        );
    }
}

/// Merge the sorted runs `from[..mid]` and `from[mid..]` into `into`.
#[allow(clippy::too_many_arguments)]
fn merge<T: Int>(
    text: &[u8],
    from: Arrays<'_, T>,
    mid: usize,
    into: Arrays<'_, T>,
    mask: Option<&SeedMask>,
    max_query_len: usize,
    n_ranges: &[Range<usize>],
    lcp_cache: &mut LcpCache,
) {
    let Arrays {
        sa: target_sa,
        lcp: target_lcp,
    } = into;
    let (mut x, mut y) = from.sa.split_at_mut(mid);
    let (mut lcp_x, mut lcp_y) = from.lcp.split_at_mut(mid);
    let mut len_x = x.len();
    let mut len_y = y.len();
    // the last LCP taken from x
    let mut m = T::default();
    let mut idx_x = 0;
    let mut idx_y = 0;
    let mut idx_target = 0;

    while idx_x < len_x && idx_y < len_y {
        let l_x = lcp_x[idx_x];

        match l_x.cmp(&m) {
            Ordering::Greater => {
                target_sa[idx_target] = x[idx_x];
                target_lcp[idx_target] = l_x;
            }
            Ordering::Less => {
                target_sa[idx_target] = y[idx_y];
                target_lcp[idx_target] = m;
                m = l_x;
            }
            Ordering::Equal => {
                let shorter_suffix = max(x[idx_x], y[idx_y]);
                let max_n = T::from_usize(text.len()) - shorter_suffix;

                let context = match mask {
                    Some(seed_mask) => T::from_usize(
                        seed_mask
                            .positions
                            .iter()
                            .filter(|&i| *i < max_n.to_usize())
                            .count(),
                    ),
                    None if max_query_len > 0 => {
                        min(T::from_usize(max_query_len), max_n)
                    }
                    None => max_n,
                };

                // LCP(X_i, Y_j)
                let (len_lcp, full_len_lcp) = if m < context {
                    let lcp = find_lcp(
                        text,
                        x[idx_x].to_usize(),
                        y[idx_y].to_usize(),
                        context - m,
                        m.to_usize(), // skip
                        mask,
                        max_query_len,
                        n_ranges,
                        Some(lcp_cache),
                    );
                    let full_lcp = match mask {
                        Some(seed_mask) => mask_offset(lcp.to_usize(), seed_mask),
                        None => lcp.to_usize(),
                    };
                    (lcp, T::from_usize(full_lcp))
                } else {
                    (context, context)
                };

                // an LCP reaching the context sorts the shorter suffix
                // first, and otherwise the byte after the LCP decides
                if len_lcp >= context {
                    target_sa[idx_target] = shorter_suffix;
                } else {
                    let cmp = text[(x[idx_x] + full_len_lcp).to_usize()]
                        .cmp(&text[(y[idx_y] + full_len_lcp).to_usize()]);

                    match cmp {
                        Ordering::Equal => {
                            target_sa[idx_target] = shorter_suffix;
                        }
                        Ordering::Less => {
                            target_sa[idx_target] = x[idx_x];
                        }
                        Ordering::Greater => {
                            target_sa[idx_target] = y[idx_y];
                        }
                    }
                }

                if target_sa[idx_target] == x[idx_x] {
                    target_lcp[idx_target] = l_x;
                } else {
                    target_lcp[idx_target] = m
                }

                m = len_lcp;
            }
        }

        if target_sa[idx_target] == x[idx_x] {
            idx_x += 1;
        } else {
            idx_y += 1;
            mem::swap(&mut x, &mut y);
            mem::swap(&mut len_x, &mut len_y);
            mem::swap(&mut lcp_x, &mut lcp_y);
            mem::swap(&mut idx_x, &mut idx_y);
        }
        idx_target += 1;
    }

    while idx_x < len_x {
        target_sa[idx_target] = x[idx_x];
        target_lcp[idx_target] = lcp_x[idx_x];
        idx_x += 1;
        idx_target += 1;
    }

    if idx_y < len_y {
        target_sa[idx_target] = y[idx_y];
        target_lcp[idx_target] = m;
        idx_y += 1;
        idx_target += 1;

        while idx_y < len_y {
            target_sa[idx_target] = y[idx_y];
            target_lcp[idx_target] = lcp_y[idx_y];
            idx_y += 1;
            idx_target += 1;
        }
    }
}

/// The LCP of each partition's first suffix with the last suffix of the
/// partition before it, or `None` for the first partition.
pub(super) fn first_lcps<T: Int>(
    text: &[u8],
    partitions: &[Partition<T>],
    mask: Option<&SeedMask>,
    max_query_len: usize,
    n_ranges: &[Range<usize>],
) -> Vec<Option<T>> {
    (0..partitions.len())
        .into_par_iter()
        .map(|i| {
            (i > 0).then(|| {
                find_lcp(
                    text,
                    partitions[i - 1].last_suffix,
                    partitions[i].first_suffix,
                    T::from_usize(text.len()),
                    0,
                    mask,
                    max_query_len,
                    n_ranges,
                    None,
                )
            })
        })
        .collect()
}

/// The LCP of the suffixes at `start1` and `start2`, which share their
/// first `skip` characters, comparing at most `len` characters.
#[allow(clippy::too_many_arguments)]
#[inline(always)]
fn find_lcp<T: Int>(
    text: &[u8],
    start1: usize,
    start2: usize,
    len: T,
    skip: usize,
    mask: Option<&SeedMask>,
    max_query_len: usize,
    n_ranges: &[Range<usize>],
    lcp_cache: Option<&mut LcpCache>,
) -> T {
    // with a mask, `len` is not consulted: the comparison
    // stops at the end of the mask or of the text
    if let Some(mask) = mask {
        return T::from_usize(
            skip + masked_lcp(
                &text[start1..],
                &text[start2..],
                &mask.positions[skip..],
            ),
        );
    }

    if let (Some(end1), Some(end2)) =
        (find_n_run(n_ranges, start1), find_n_run(n_ranges, start2))
    {
        let run = min(end1 - start1, end2 - start2);
        return T::from_usize(if max_query_len > 0 {
            min(run, max_query_len)
        } else {
            run
        });
    }

    // `len` counts bytes past `skip`, so the cap leaves
    // max_query_len - skip of them
    let len = if max_query_len > 0 {
        max_query_len - skip
    } else {
        len.to_usize()
    };
    T::from_usize(suffix_lcp(text, start1, start2, skip, len, lcp_cache))
}

/// The end of the run of `N`s in `n_ranges` holding `suffix`, if it is
/// in one.
fn find_n_run(n_ranges: &[Range<usize>], suffix: usize) -> Option<usize> {
    n_ranges
        .binary_search_by(|range| {
            if range.contains(&suffix) {
                Ordering::Equal
            } else if range.start < suffix {
                Ordering::Less
            } else {
                Ordering::Greater
            }
        })
        .ok()
        .map(|i| n_ranges[i].end)
}

#[cfg(test)]
mod test {
    use super::find_lcp;
    use crate::format::SeedMask;
    use anyhow::Result;
    use pretty_assertions::assert_eq;

    #[test]
    fn test_find_lcp_no_seed_mask() -> Result<()> {
        //           012345
        let text = b"TTTAGC";

        // 0: TTTAGC
        // 1:  TTAGC
        // 6: len of text
        assert_eq!(find_lcp::<u32>(text, 0, 1, 6, 0, None, 0, &[], None), 2);

        // 0: TTTAGC
        // 2:   TAGC
        // 6: len of text
        assert_eq!(find_lcp::<u32>(text, 0, 2, 6, 0, None, 0, &[], None), 1);

        // 0: TTTAGC
        // 1:  TTAGC
        // 1: max query len = 1
        assert_eq!(find_lcp::<u32>(text, 0, 1, 1, 0, None, 0, &[], None), 1);

        // 0: TTTAGC
        // 3:    AGC
        // 6: len of text
        assert_eq!(find_lcp::<u32>(text, 0, 3, 6, 0, None, 0, &[], None), 0);

        // 0: TTTTAG
        // 1:  TTTAG
        // skip 2 with max query len 2: nothing left to compare
        let text = b"TTTTAG";
        assert_eq!(find_lcp::<u32>(text, 0, 1, 6, 2, None, 2, &[], None), 2);

        // skip 1 with max query len 3: one byte left, T == T
        assert_eq!(find_lcp::<u32>(text, 0, 1, 6, 1, None, 3, &[], None), 3);

        // both suffixes inside a run of four N, max query len 2
        let text = b"NNNNA";
        assert_eq!(find_lcp::<u32>(text, 0, 1, 5, 0, None, 2, &[0..4], None), 2);

        Ok(())
    }

    #[test]
    fn test_find_lcp_with_seed_mask() -> Result<()> {
        //           012345
        let text = b"TTTTTA";
        let mask = SeedMask::new("1101")?;
        let mask = Some(&mask);

        // 0: TTTTTA
        // 1:  TTTTA
        assert_eq!(find_lcp::<u32>(text, 0, 1, 3, 0, mask, 0, &[], None), 3);

        // 0: TTTTTA
        // 2:   TTTA
        assert_eq!(find_lcp::<u32>(text, 0, 2, 3, 0, mask, 0, &[], None), 2);

        // 0: TTTTTA
        // 5:      A
        assert_eq!(find_lcp::<u32>(text, 0, 5, 3, 0, mask, 0, &[], None), 0);

        Ok(())
    }
}
