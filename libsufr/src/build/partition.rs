use super::scratch::{Scratch, ScratchKind};
use crate::{format::SeedMask, io::Int};
use anyhow::{bail, Result};
use log::info;
use rayon::prelude::*;
use std::{sync::Mutex, time::Instant};

const NUCLEIC_ACID: &[u8] = b"$%ACGNT";
const AMINO_ACID: &[u8] = b"$%*-ABCDEFGHIJKLMNOPQRSTUVWXYZ";

/// The characters a suffix is placed by, ranked in sort order, and how
/// many of them a bucket number packs into a `u16`.
pub(super) struct PrefixEncoder {
    /// The rank of each byte.
    table: [u16; 256],

    /// The number of bits per character in a bucket number.
    bits: usize,

    /// The number of leading characters in a bucket number.
    count: usize,
}

impl PrefixEncoder {
    fn nucleic_acid() -> Self {
        Self::ranked(NUCLEIC_ACID)
    }

    fn amino_acid() -> Self {
        Self::ranked(AMINO_ACID)
    }

    fn bytes() -> Self {
        let mut table = [0u16; 256];
        for (c, rank) in table.iter_mut().enumerate() {
            *rank = c as u16;
        }
        Self {
            table,
            bits: 8,
            count: 2,
        }
    }

    /// The encoder for `characters`, which are in ascending order.
    fn ranked(characters: &[u8]) -> Self {
        let mut table = [0u16; 256];
        for (i, &c) in characters.iter().enumerate() {
            table[c as usize] = i as u16;
        }
        let bits = characters.len().next_power_of_two().trailing_zeros() as usize;
        Self {
            table,
            bits,
            count: u16::BITS as usize / bits,
        }
    }

    /// The nucleic acid encoder when `occupancy` marks only its bytes,
    /// else the amino acid encoder when it marks only those, else bytes.
    pub(super) fn from_occupancy(mut occupancy: [bool; 256]) -> Self {
        for &c in NUCLEIC_ACID {
            occupancy[c as usize] = false;
        }
        if occupancy.iter().all(|&seen| !seen) {
            info!("Detected Nucleic Acid alphabet: using 5-character prefix for partitioning");
            return Self::nucleic_acid();
        }

        for &c in AMINO_ACID {
            occupancy[c as usize] = false;
        }
        if occupancy.iter().all(|&seen| !seen) {
            info!("Detected Amino Acid alphabet: using 3-character prefix for partitioning");
            Self::amino_acid()
        } else {
            info!("Detected non-Nucleic/non-Amino Acid: using 2-byte prefix for partitioning");
            Self::bytes()
        }
    }

    /// The number of distinct bucket numbers.
    pub(super) fn num_buckets(&self) -> usize {
        1 << (self.bits * self.count)
    }

    /// The bucket number of the suffix at `i`: its first `count`
    /// characters, or those at the care positions of `mask`.
    #[inline(always)]
    pub(super) fn encode(&self, text: &[u8], i: usize, mask: Option<&SeedMask>) -> u16 {
        let (bits, count) = (self.bits, self.count);
        let mut key = 0;
        match mask {
            Some(mask) => {
                for (j, &off) in mask.positions.iter().enumerate().take(count) {
                    let c = text.get(i + off).copied().unwrap_or(0);
                    key |= self.table[c as usize] << ((count - j - 1) * bits);
                }
            }
            None => {
                for (j, &c) in text[i..].iter().enumerate().take(count) {
                    key |= self.table[c as usize] << ((count - j - 1) * bits);
                }
            }
        }
        key
    }
}

/// The sorted suffixes and LCPs of one partition.
#[derive(Debug)]
pub(super) struct Partition<T: Int> {
    /// The position of this partition in the suffix array.
    pub(super) order: usize,

    /// The number of suffixes/LCP values contained in this partition.
    pub(super) len: usize,

    /// The first suffix, for the LCP at the boundary with the partition
    /// before.
    pub(super) first_suffix: usize,

    /// The last suffix, for the LCP at the boundary with the partition
    /// after.
    pub(super) last_suffix: usize,

    /// The buffer containing the suffix array.
    pub(super) sa_buffer: Scratch<T>,

    /// The buffer containing the LCP array, when one is written.
    pub(super) lcp_buffer: Option<Scratch<T>>,
}

/// The unsorted suffixes in one scratch buffer per prefix, in prefix
/// order.
#[derive(Debug)]
pub(super) struct PartitionedSuffixes<T: Int> {
    buffers: Vec<Scratch<T>>,
    pub(super) num_suffixes: usize,
}

/// Place each suffix of `text` whose first byte `keep` marks in the
/// scratch buffer of `kind` for its bucket number from `encoder`.
pub(super) fn partition<T: Int>(
    text: &[u8],
    encoder: &PrefixEncoder,
    keep: &[bool; 256],
    mask: Option<&SeedMask>,
    kind: ScratchKind,
) -> Result<PartitionedSuffixes<T>> {
    let num_partitions = encoder.num_buckets();
    let mut buffers: Vec<_> = vec![];
    for _ in 0..num_partitions {
        buffers.push(Mutex::new(Scratch::new(kind)));
    }

    let now = Instant::now();
    const CHUNK_SIZE: usize = 1 << 19;
    (0..text.len().div_ceil(CHUNK_SIZE))
        .into_par_iter()
        .try_for_each(|chunk| -> Result<()> {
            let start = chunk * CHUNK_SIZE;
            let end = (start + CHUNK_SIZE).min(text.len());

            let mut suffixes = Vec::with_capacity(end - start);
            let mut keys = Vec::with_capacity(end - start);
            let mut counts = vec![0; num_partitions];

            for i in start..end {
                if !keep[text[i] as usize] {
                    continue;
                }

                let key = encoder.encode(text, i, mask);
                keys.push(key);
                counts[key as usize] += 1;
                suffixes.push(T::from_usize(i));
            }

            let mut offsets = Vec::with_capacity(num_partitions);
            let mut total = 0;
            for &n in counts.iter() {
                offsets.push(total);
                total += n;
            }

            let mut ordered = vec![T::default(); keys.len()];
            let mut cursors = offsets.clone();
            for (&key, &suf) in keys.iter().zip(suffixes.iter()) {
                let offset = &mut cursors[key as usize];
                ordered[*offset] = suf;
                *offset += 1;
            }

            for (partition_num, (&offset, &count)) in
                offsets.iter().zip(counts.iter()).enumerate()
            {
                if count > 0 {
                    let slice = &ordered[offset..offset + count];
                    match buffers[partition_num].lock() {
                        Ok(mut buf) => {
                            if buf.extend_from_slice(slice).is_err() {
                                bail!("Unable to write data to disk")
                            }
                        }
                        Err(e) => bail!("{e}"),
                    }
                }
            }
            Ok(())
        })?;

    let mut num_suffixes = 0;
    let buffers = buffers
        .into_iter()
        .map(|buffer| match buffer.into_inner() {
            Ok(mut buf) => {
                buf.flush()?;
                num_suffixes += buf.count();
                Ok(buf)
            }
            Err(e) => panic!("Failed to lock: {e}"),
        })
        .collect::<Result<Vec<_>>>()?;

    info!(
        "Wrote {num_suffixes} unsorted suffixes to partition{} in {:?}",
        if num_partitions == 1 { "" } else { "s" },
        now.elapsed()
    );

    Ok(PartitionedSuffixes {
        buffers,
        num_suffixes,
    })
}

/// Gather the buffers, which are in prefix order, into `num_partitions`
/// groups of about the same number of suffixes.
pub(super) fn group_buffers<T: Int>(
    partitioned: PartitionedSuffixes<T>,
    num_partitions: usize,
) -> Result<Vec<Vec<Scratch<T>>>> {
    // round up to take every suffix
    let num_per_partition =
        (partitioned.num_suffixes as f64 / num_partitions as f64).ceil() as usize;
    let mut num_taken = 0;
    let mut partition_inputs = std::iter::repeat_with(Vec::new)
        .take(num_partitions)
        .collect::<Vec<_>>();

    // there are (probably) many more buffers than partitions,
    // so each partition takes buffers from the left until it
    // passes its share, which evens out the sorting work
    let mut part_buffers = partitioned.buffers.into_iter();
    for (partition_num, partition_input) in partition_inputs.iter_mut().enumerate() {
        let boundary = num_per_partition * (partition_num + 1);
        for buffer in part_buffers.by_ref() {
            if buffer.count() > 0 {
                num_taken += buffer.count();
                partition_input.push(buffer);
            }

            // the last partition takes the rest
            if partition_num < num_partitions - 1 && num_taken > boundary {
                break;
            }
        }
    }

    if num_taken != partitioned.num_suffixes {
        bail!(
            "Took {num_taken} but needed to take {}",
            partitioned.num_suffixes
        );
    }
    Ok(partition_inputs)
}

#[cfg(test)]
mod test {
    use super::PrefixEncoder;
    use crate::format::SeedMask;
    use anyhow::Result;

    #[test]
    fn test_encoder_properties() {
        let nucleic_acid = PrefixEncoder::nucleic_acid();
        assert_eq!(nucleic_acid.bits, 3);
        assert_eq!(nucleic_acid.count, 5);

        let amino_acid = PrefixEncoder::amino_acid();
        assert_eq!(amino_acid.bits, 5);
        assert_eq!(amino_acid.count, 3);
    }

    #[test]
    fn test_encode() -> Result<()> {
        // ranks in $%ACGNT: $=0 %=1 A=2 C=3 G=4 N=5 T=6
        let text = b"ACGTNA$";
        let mut occupancy = [false; 256];
        for &b in text {
            occupancy[b as usize] = true;
        }
        let encoder = PrefixEncoder::from_occupancy(occupancy);
        assert_eq!(encoder.num_buckets(), 1 << 15);
        assert_eq!(encoder.encode(text, 0, None), 0o23465);

        // a suffix shorter than 5 leaves the low places 0
        assert_eq!(encoder.encode(text, 4, None), 0o52000);

        let mask = SeedMask::new("1101")?;
        assert_eq!(encoder.encode(text, 0, Some(&mask)), 0o23600);
        assert_eq!(encoder.encode(text, 5, Some(&mask)), 0o20000);
        Ok(())
    }
}
