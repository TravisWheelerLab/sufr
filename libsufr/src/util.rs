//! Utility functions

use crate::types::{
    FromUsize, Int, SequenceFileData, SuffixSortType, OUTFILE_VERSION,
    SENTINEL_CHARACTER,
};
use anyhow::{anyhow, bail, Result};
use needletail::parse_fastx_file;
use std::{fs::File, io::Read, path::Path, slice};

// --------------------------------------------------
/// When using a seed mask, the LCP stored on disk is the number of "care"
/// positions held in common, but the *actual* LCP between two strings
/// should include any "don't-care" positions up to the next "care" position.
///
/// Args:
/// * `lcp`: the LCP (probably from LCP array on disk)
/// * `sort_type`: the suffixes are fully sorted or masked
pub fn find_lcp_full_offset(lcp: usize, sort_type: &SuffixSortType) -> usize {
    match sort_type {
        SuffixSortType::Mask(seed_mask) => {
            if lcp == 0 || lcp > seed_mask.bytes.len() {
                lcp
            } else {
                // E.g., LCP = 1, so get the 0th offset
                let offset = seed_mask.positions[lcp - 1];
                let next_offset = *seed_mask.positions.get(lcp).unwrap_or(&0);
                if (next_offset > offset) && (next_offset - offset) > 1 {
                    next_offset
                } else {
                    offset + 1
                }
            }
        }
        _ => lcp,
    }
}

// --------------------------------------------------
/// Read sequence data from file
///
/// Args:
/// * `path`: to FASTA/FASTQ file
/// * `sequence_delimiter`: the character/byte you wish to place
///    between sequences. For many applications (both nucleotide and
///    protein), you might choose a character like `%` that sorts below
///    the alphabet (A-Z) but above the sentinel `$` that is placed at
///    the end of the returned text; however, when building a suffix array 
///    for use in a Burrows-Wheeler Transform, it might be better to use 
///    `N` for nucleotides and `X` for protein.
pub fn read_sequence_file(
    path: &Path,
    sequence_delimiter: u8,
) -> Result<SequenceFileData> {
    let mut reader = parse_fastx_file(path)?;
    let mut seq: Vec<u8> = Vec::with_capacity(u32::MAX as usize);
    let mut sequence_names: Vec<String> = vec![];
    let mut start_positions: Vec<usize> = vec![];
    let mut i = 0;
    while let Some(rec) = reader.next() {
        let rec = rec?;
        if i > 0 {
            seq.push(sequence_delimiter);
        }

        // Record current length as start position
        start_positions.push(seq.len());
        let mut tmp: Vec<u8> = rec.seq().iter().copied().collect();
        seq.append(&mut tmp);
        i += 1;

        // Only take ID value up to first whitespace
        let id = String::from_utf8(rec.id().to_vec())?
            .split_whitespace()
            .next()
            .map_or((i + 1).to_string(), |v| v.to_string());

        sequence_names.push(id);
    }

    // File delimiter
    seq.push(SENTINEL_CHARACTER);

    Ok(SequenceFileData {
        seq,
        start_positions,
        sequence_names,
    })
}

// --------------------------------------------------
/// Find length of the input text from a _.sufr_
/// file to determine the `Int` type, `u32` or `u64`
///
/// Args:
/// * `filename`: the _.sufr_ filename to read
pub fn read_text_length(filename: &str) -> Result<usize> {
    let mut file = File::open(filename).map_err(|e| anyhow!("{filename}: {e}"))?;

    // Meta (version, is_dna)
    let mut buffer = [0; 4];
    file.read_exact(&mut buffer)?;

    let outfile_version = buffer[0];
    if outfile_version == OUTFILE_VERSION {
        // Length of text is the next usize
        let mut buffer = [0; 8];
        file.read_exact(&mut buffer)?;
        Ok(usize::from_ne_bytes(buffer))
    } else {
        bail!("Unknown sufr version {outfile_version}");
    }
}

// --------------------------------------------------
/// Convert a slice of raw U8 read from disk into a
/// `Vec<T>` (where `T` is the `Int` 32/64)
///
/// Args:
/// * `buffer`: vector of raw `u8` values (from disk)
/// * `len`: the number of `T` values in the resulting vector
pub fn slice_u8_to_vec<T>(buffer: &[u8], len: usize) -> Vec<T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + serde::ser::Serialize,
{
    unsafe { std::slice::from_raw_parts(buffer.as_ptr() as *const _, len).to_vec() }
}

// --------------------------------------------------
/// Turn a `usize` into a vector of `u8` for serializing to disk.
///
/// Args:
/// * `value`: a `usize`
pub fn usize_to_bytes(value: usize) -> Vec<u8> {
    // Determine the size of usize in bytes
    let size = std::mem::size_of::<usize>();

    // Create a vector to hold the bytes
    let mut bytes = Vec::with_capacity(size);

    // Convert usize to bytes
    for i in 0..size {
        bytes.push((value >> (i * 8)) as u8);
    }

    bytes
}

// --------------------------------------------------
/// Convert a `Vec<T>` (where `T` is the `Int` 32/64) into a
/// slice of raw U8 for serializing to disk
///
/// Args:
/// * `vec`: a vector of `T` values
pub fn vec_to_slice_u8<T>(vec: &[T]) -> &[u8]
where
    T: Int + FromUsize<T> + Sized + Send + Sync + serde::ser::Serialize,
{
    unsafe {
        slice::from_raw_parts(
            vec.as_ptr() as *const _,
            std::mem::size_of_val(vec), //vec.len() * mem::size_of::<T>(),
        )
    }
}

// --------------------------------------------------
#[cfg(test)]
mod tests {
    use super::{
        find_lcp_full_offset, read_sequence_file, read_text_length, slice_u8_to_vec,
        usize_to_bytes, vec_to_slice_u8,
    };
    use crate::types::{SeedMask, SuffixSortType};
    use anyhow::Result;
    use pretty_assertions::assert_eq;
    use std::path::Path;

    #[test]
    fn test_read_sequence_file() -> Result<()> {
        let file = Path::new("../data/inputs/2.fa");
        let sequence_delimiter = b'N';
        let res = read_sequence_file(&file, sequence_delimiter);
        assert!(res.is_ok());
        let data = res.unwrap();
        assert_eq!(data.seq, b"ACGTacgtNacgtACGT$");
        assert_eq!(data.start_positions, [0, 9]);
        assert_eq!(data.sequence_names, ["ABC", "DEF"]);
        Ok(())
    }

    #[test]
    fn test_read_text_length() -> Result<()> {
        let sufr_file = "../data/inputs/2.sufr";
        let res = read_text_length(sufr_file);
        assert!(res.is_ok());
        let len = res.unwrap();
        assert_eq!(len, 18);
        Ok(())
    }

    #[test]
    fn test_usize_to_bytes() -> Result<()> {
        assert_eq!(usize_to_bytes(usize::MIN), [0, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(usize_to_bytes(1), [1, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(usize_to_bytes(10), [10, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(usize_to_bytes(100), [100, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(usize_to_bytes(1000), [232, 3, 0, 0, 0, 0, 0, 0]);
        assert_eq!(usize_to_bytes(10000), [16, 39, 0, 0, 0, 0, 0, 0]);
        assert_eq!(usize_to_bytes(100000), [160, 134, 1, 0, 0, 0, 0, 0]);
        assert_eq!(
            usize_to_bytes(usize::MAX),
            [255, 255, 255, 255, 255, 255, 255, 255]
        );
        Ok(())
    }

    #[test]
    fn test_slice_u8_to_vec() -> Result<()> {
        let res: Vec<u32> = slice_u8_to_vec(&[0, 0, 0, 0], 1);
        assert_eq!(res, &[0u32]);

        let res: Vec<u64> = slice_u8_to_vec(&[0, 0, 0, 0, 0, 0, 0, 0], 1);
        assert_eq!(res, &[0u64]);

        let res: Vec<u32> = slice_u8_to_vec(&[1, 0, 0, 0], 1);
        assert_eq!(res, &[1u32]);

        let res: Vec<u64> = slice_u8_to_vec(&[1, 0, 0, 0, 0, 0, 0, 0], 1);
        assert_eq!(res, &[1u64]);

        let res: Vec<u32> = slice_u8_to_vec(&[255, 255, 255, 255], 1);
        assert_eq!(res, &[u32::MAX]);

        let res: Vec<u64> =
            slice_u8_to_vec(&[255, 255, 255, 255, 255, 255, 255, 255], 1);
        assert_eq!(res, &[u64::MAX]);

        let res: Vec<u32> = slice_u8_to_vec(&[0, 0, 0, 0, 255, 255, 255, 255], 2);
        assert_eq!(res, &[0u32, u32::MAX]);

        let res: Vec<u64> = slice_u8_to_vec(
            &[
                0, 0, 0, 0, 0, 0, 0, 0, 255, 255, 255, 255, 255, 255, 255, 255,
            ],
            2,
        );
        assert_eq!(res, &[0u64, u64::MAX]);

        Ok(())
    }

    #[test]
    fn test_vec_to_slice_u8() -> Result<()> {
        let res = vec_to_slice_u8(&[0u32]);
        assert_eq!(res, &[0, 0, 0, 0]);

        let res = vec_to_slice_u8(&[0u64]);
        assert_eq!(res, &[0, 0, 0, 0, 0, 0, 0, 0]);

        let res = vec_to_slice_u8(&[1u32]);
        assert_eq!(res, &[1, 0, 0, 0]);

        let res = vec_to_slice_u8(&[1u64]);
        assert_eq!(res, &[1, 0, 0, 0, 0, 0, 0, 0]);

        let res = vec_to_slice_u8(&[u32::MAX]);
        assert_eq!(res, &[255, 255, 255, 255]);

        let res = vec_to_slice_u8(&[u64::MAX]);
        assert_eq!(res, &[255, 255, 255, 255, 255, 255, 255, 255]);

        let res = vec_to_slice_u8(&[0u32, u32::MAX]);
        assert_eq!(res, &[0, 0, 0, 0, 255, 255, 255, 255]);

        let res = vec_to_slice_u8(&[0u64, u64::MAX]);
        assert_eq!(
            res,
            &[0, 0, 0, 0, 0, 0, 0, 0, 255, 255, 255, 255, 255, 255, 255, 255]
        );

        Ok(())
    }

    #[test]
    fn test_find_lcp_full_offset() -> Result<()> {
        let seed_mask = SeedMask::new("101")?;
        let sort_type = SuffixSortType::Mask(seed_mask);
        assert_eq!(find_lcp_full_offset(0, &sort_type), 0);
        assert_eq!(find_lcp_full_offset(1, &sort_type), 2);
        assert_eq!(find_lcp_full_offset(2, &sort_type), 3);

        let seed_mask = SeedMask::new("11011")?;
        let sort_type = SuffixSortType::Mask(seed_mask);
        assert_eq!(find_lcp_full_offset(0, &sort_type), 0);
        assert_eq!(find_lcp_full_offset(1, &sort_type), 1);
        assert_eq!(find_lcp_full_offset(2, &sort_type), 3);
        assert_eq!(find_lcp_full_offset(3, &sort_type), 4);
        assert_eq!(find_lcp_full_offset(4, &sort_type), 5);

        let seed_mask = SeedMask::new("10011001")?;
        let sort_type = SuffixSortType::Mask(seed_mask);
        assert_eq!(find_lcp_full_offset(0, &sort_type), 0);
        assert_eq!(find_lcp_full_offset(1, &sort_type), 3);
        assert_eq!(find_lcp_full_offset(2, &sort_type), 4);
        assert_eq!(find_lcp_full_offset(3, &sort_type), 7);
        assert_eq!(find_lcp_full_offset(4, &sort_type), 8);

        Ok(())
    }
}
