//! Utility functions

use crate::types::{SequenceFileData, SuffixSortType, SENTINEL_CHARACTER};
use anyhow::Result;
use needletail::parse_fastx_file;
use std::path::Path;

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
///   between sequences. For many applications (both nucleotide and
///   protein), you might choose a character like `%` that sorts below
///   the alphabet (A-Z) but above the sentinel `$` that is placed at
///   the end of the returned text; however, when building a suffix array
///   for use in a Burrows-Wheeler Transform, it might be better to use
///   `N` for nucleotides and `X` for protein.
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
        seq.extend_from_slice(&rec.seq());
        i += 1;

        // Only take ID value up to first whitespace
        let id = std::str::from_utf8(rec.id())?
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
#[cfg(test)]
mod tests {
    use super::{find_lcp_full_offset, read_sequence_file};
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
