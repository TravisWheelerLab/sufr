use super::partition::Partition;
use crate::{
    format::Header,
    io::{as_bytes, write_all_at, Int},
};
use anyhow::{anyhow, Result};
use log::info;
use rayon::prelude::*;
use std::{fs::File, mem, path::Path, time::Instant};

/// Write `text`, then each partition's suffix array (and LCP array, when
/// `header` places one, starting with its entry in `first_lcps`) in
/// parallel at its offset, then `header` and the sequence names.
pub(super) fn write<T: Int>(
    path: &Path,
    header: &Header,
    text: &[u8],
    partitions: Vec<Partition<T>>,
    first_lcps: &[Option<T>],
) -> Result<()> {
    let now = Instant::now();
    let int_size = mem::size_of::<T>();
    let file = File::create(path).map_err(|e| anyhow!("{}: {e}", path.display()))?;
    write_all_at(&file, text, header.text_pos as u64)?;

    let offsets = partition_offsets(&partitions);

    // consuming the partitions' buffers frees memory
    // and disk space along the way
    partitions
        .into_par_iter()
        .zip(offsets)
        .enumerate()
        .try_for_each(|(i, (partition, offset))| -> Result<()> {
            let sa = partition.sa_buffer.consume()?;
            write_all_at(
                &file,
                as_bytes(&sa),
                (header.suffix_array_pos + offset * int_size) as u64,
            )?;
            drop(sa);

            if let (Some(lcp_buffer), Some(lcp_pos)) =
                (partition.lcp_buffer, header.lcp_pos)
            {
                let mut lcp = lcp_buffer.consume()?;
                if let Some(first_lcp) = first_lcps[i] {
                    if let Some(val) = lcp.first_mut() {
                        *val = first_lcp;
                    }
                }
                write_all_at(
                    &file,
                    as_bytes(&lcp),
                    (lcp_pos + offset * int_size) as u64,
                )?;
            }
            Ok(())
        })?;

    header.write(&file)?;
    info!(
        "Wrote suffix array to '{}' in {:?}",
        path.display(),
        now.elapsed()
    );
    Ok(())
}

/// The offset (in elements) of each partition in the arrays.
fn partition_offsets<T: Int>(partitions: &[Partition<T>]) -> Vec<usize> {
    let mut offsets = Vec::with_capacity(partitions.len());
    let mut total = 0;
    for partition in partitions {
        offsets.push(total);
        total += partition.len;
    }
    offsets
}
