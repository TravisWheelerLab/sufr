use super::partition::Partition;
use crate::{
    buffer::Buffer,
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
pub(super) fn write<W: Int, S: Buffer<W>>(
    path: &Path,
    header: &Header,
    text: &[u8],
    partitions: Vec<Partition<S>>,
    first_lcps: &[Option<W>],
) -> Result<()> {
    let now = Instant::now();
    let int_size = mem::size_of::<W>();
    let file = File::create(path).map_err(|e| anyhow!("{}: {e}", path.display()))?;
    write_all_at(&file, text, header.text_pos as u64)?;

    let offsets = partition_offsets(&partitions);

    // dropping each buffer once it is written frees memory
    // and disk space along the way
    partitions
        .into_par_iter()
        .zip(offsets)
        .enumerate()
        .try_for_each(|(i, (partition, offset))| -> Result<()> {
            write_all_at(
                &file,
                as_bytes(&partition.sa_buffer.range(0..partition.len)?),
                (header.suffix_array_pos + offset * int_size) as u64,
            )?;
            drop(partition.sa_buffer);

            if let (Some(lcp_buffer), Some(lcp_pos)) =
                (partition.lcp_buffer, header.lcp_pos)
            {
                let mut lcp = lcp_buffer.range(0..partition.len)?.into_owned();
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
fn partition_offsets<W: Int, S: Buffer<W>>(partitions: &[Partition<S>]) -> Vec<usize> {
    let mut offsets = Vec::with_capacity(partitions.len());
    let mut total = 0;
    for partition in partitions {
        offsets.push(total);
        total += partition.len;
    }
    offsets
}
