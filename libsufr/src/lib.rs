// #![warn(missing_docs)]

//! # Parallel Construction of Suffix Arrays in Rust
//!
//! Inspired by the merge sort approach described in "Cache-friendly,
//! Parallel, and Samplesort-based Constructor for Suffix Arrays and LCP
//! Arrays." [^patro]
//!
//! [^patro]: <https://doi.org/10.4230/LIPIcs.WABI.2023.16>
//!
//! * To create a new suffix array/_.sufr_ file, see [sufr_builder]
//! * To interact and query an existing suffix array/_.sufr_ file, see [sufr_file]
//!
//! ## Authors
//! * Ken Youens-Clark <kyclark@gmail.com>
//! * Jack Roddy <jroddy@arizona.edu>
//! * Travis Wheeler <twheeler@arizona.edu>

mod file_access;
pub mod suffix_array;
pub mod sufr_builder;
pub mod sufr_file;
mod sufr_search;
pub mod types;
pub mod util;

// --------------------------------------------------
#[cfg(test)]
mod tests {
    use super::{
        sufr_builder::SufrBuilder,
        sufr_file::SufrFile,
        types::{SufrBuilderArgs, OUTFILE_VERSION},
        util::read_sequence_file,
    };
    use anyhow::Result;
    use std::path::Path;
    use tempfile::NamedTempFile;

    #[test]
    fn test_write_read_suffix_file_32() -> Result<()> {
        let seq_file = Path::new("../data/inputs/2.fa");
        let sequence_delimiter = b'N';
        let seq_data = read_sequence_file(seq_file, sequence_delimiter)?;
        let outfile = NamedTempFile::new()?;
        let outpath = outfile.path().to_string_lossy().to_string();
        let args = SufrBuilderArgs {
            text: seq_data.seq,
            low_memory: true,
            path: Some(outpath.clone()),
            max_query_len: None,
            is_dna: true,
            allow_ambiguity: false,
            ignore_softmask: false,
            sequence_starts: seq_data.start_positions.into_iter().collect(),
            sequence_names: seq_data.sequence_names,
            num_partitions: 2,
            seed_mask: None,
            random_seed: 0,
        };
        let res = SufrBuilder::<u32>::new(args);
        assert!(res.is_ok());

        let builder = res.unwrap();
        assert!(Path::new(&builder.path).exists());

        let res: Result<SufrFile<u32>> = SufrFile::read(&outpath, false);
        assert!(res.is_ok());

        let mut sufr_file = res.unwrap();
        assert_eq!(sufr_file.version, OUTFILE_VERSION);
        assert!(sufr_file.is_dna);
        assert_eq!(sufr_file.text_len, 18);
        assert_eq!(sufr_file.num_sequences, 2);
        assert_eq!(sufr_file.sequence_starts, [0, 9]);
        assert_eq!(sufr_file.sequence_names, ["ABC", "DEF"]);
        assert_eq!(sufr_file.text, b"ACGTACGTNACGTACGT$");

        let file_sa: Vec<_> = sufr_file.suffix_array_file.iter().collect();
        let sorted_sa = [17, 13, 9, 0, 4, 14, 10, 1, 5, 15, 11, 2, 6, 16, 12, 3, 7];
        assert_eq!(file_sa, sorted_sa);

        let file_lcp: Vec<_> = sufr_file.lcp_file.iter().collect();
        let lcp = [0, 0, 4, 8, 4, 0, 3, 7, 3, 0, 2, 6, 2, 0, 1, 5, 1];
        assert_eq!(file_lcp, lcp);
        Ok(())
    }

    #[test]
    fn test_write_read_suffix_file_64() -> Result<()> {
        let seq_file = Path::new("../data/inputs/1.fa");
        let sequence_delimiter = b'N';
        let seq_data = read_sequence_file(seq_file, sequence_delimiter)?;
        let outfile = NamedTempFile::new()?;
        let outpath = outfile.path().to_string_lossy().to_string();
        let args = SufrBuilderArgs {
            text: seq_data.seq,
            low_memory: true,
            path: Some(outpath.clone()),
            max_query_len: None,
            is_dna: true,
            allow_ambiguity: true,
            ignore_softmask: false,
            sequence_starts: seq_data.start_positions.into_iter().collect(),
            sequence_names: seq_data.sequence_names,
            num_partitions: 2,
            seed_mask: None,
            random_seed: 0,
        };

        let res = SufrBuilder::<u64>::new(args);
        assert!(res.is_ok());

        let builder = res.unwrap();
        assert!(Path::new(&builder.path).exists());

        let res: Result<SufrFile<u64>> = SufrFile::read(&outpath, false);
        assert!(res.is_ok());

        let mut sufr_file = res.unwrap();
        assert_eq!(sufr_file.version, OUTFILE_VERSION);
        assert!(sufr_file.is_dna);
        assert_eq!(sufr_file.text_len, 11);
        assert_eq!(sufr_file.num_sequences, 1);
        assert_eq!(sufr_file.sequence_starts, [0]);
        assert_eq!(sufr_file.sequence_names, ["1"]);
        assert_eq!(sufr_file.text, b"ACGTNNACGT$");
        assert_eq!(sufr_file.len_suffixes, 11);

        let file_sa: Vec<_> = sufr_file.suffix_array_file.iter().collect();
        assert_eq!(file_sa, &[10, 6, 0, 7, 1, 8, 2, 5, 4, 9, 3]);
        let file_lcp: Vec<_> = sufr_file.lcp_file.iter().collect();
        assert_eq!(file_lcp, &[0, 0, 4, 0, 3, 0, 2, 0, 1, 0, 1]);
        Ok(())
    }

    #[test]
    fn test_subsample_suffix_array() -> Result<()> {
        let seq_file = Path::new("../data/inputs/smol.fa");
        let sequence_delimiter = b'N';
        let seq_data = read_sequence_file(seq_file, sequence_delimiter)?;
        let outfile = NamedTempFile::new()?;
        let outpath = outfile.path().to_string_lossy().to_string();
        let builder_args = SufrBuilderArgs {
            text: seq_data.seq,
            low_memory: true,
            path: Some(outpath.clone()),
            max_query_len: None,
            is_dna: true,
            allow_ambiguity: false,
            ignore_softmask: false,
            sequence_starts: seq_data.start_positions,
            sequence_names: seq_data.sequence_names,
            num_partitions: 2,
            seed_mask: None,
            random_seed: 0,
        };
        let res = SufrBuilder::<u32>::new(builder_args);
        assert!(res.is_ok());

        let builder = res.unwrap();
        assert_eq!(builder.num_suffixes, 364);

        let mut sufr_file: SufrFile<u32> = SufrFile::read(&outpath, false)?;
        let full_sa: Vec<_> = sufr_file.suffix_array_file.iter().collect();
        let full_lcp: Vec<_> = sufr_file.lcp_file.iter().collect();
        assert_eq!(full_sa.len(), 364);
        assert_eq!(full_lcp.len(), 364);

        let max_query_len = 1;
        let (sub_sa, sub_rank) = sufr_file.subsample_suffix_array(max_query_len);
        assert_eq!(sub_sa.len(), 5);
        assert_eq!(sub_rank.len(), 5);
        // $, A, C, G, T
        assert_eq!(sub_sa, vec![365, 364, 92, 224, 363]);
        assert_eq!(sub_rank, vec![0, 1, 94, 191, 284]);

        let max_query_len = 2;
        let (sub_sa, sub_rank) = sufr_file.subsample_suffix_array(max_query_len);
        assert_eq!(sub_sa.len(), 20);
        assert_eq!(sub_rank.len(), 20);
        // $, A$, AA, AC, AG, AN, AT, CA, CC, CG, CT, GA, GC, GG, GN, GT, TA, TC, TG, TT
        assert_eq!(
            sub_sa,
            vec![
                365, 364, 358, 91, 341, 255, 362, 92, 339, 233, 296, 224, 88, 129, 110,
                96, 363, 217, 223, 356
            ]
        );
        assert_eq!(
            sub_rank,
            vec![
                0, 1, 2, 38, 49, 70, 71, 94, 112, 143, 170, 191, 216, 252, 269, 270,
                284, 298, 315, 343
            ]
        );

        let max_query_len = 3;
        let (sub_sa, sub_rank) = sufr_file.subsample_suffix_array(max_query_len);
        assert_eq!(sub_sa.len(), 71);
        assert_eq!(sub_rank.len(), 71);

        let max_query_len = 5;
        let (sub_sa, sub_rank) = sufr_file.subsample_suffix_array(max_query_len);
        assert_eq!(sub_sa.len(), 293);
        assert_eq!(sub_rank.len(), 293);

        let max_query_len = 3;
        let (sub_sa, sub_rank) = sufr_file.subsample_suffix_array(max_query_len);
        assert_eq!(sub_sa.len(), 71);
        assert_eq!(sub_rank.len(), 71);

        Ok(())
    }

    #[test]
    fn test_spaced_seeds_1() -> Result<()> {
        let seq_file = Path::new("../data/inputs/mostlya1.fa");
        let sequence_delimiter = b'N';
        let seq_data = read_sequence_file(seq_file, sequence_delimiter)?;
        let outfile = NamedTempFile::new()?;
        let outpath = outfile.path().to_string_lossy().to_string();
        let builder_args = SufrBuilderArgs {
            text: seq_data.seq,
            low_memory: true,
            path: Some(outpath.clone()),
            max_query_len: None,
            is_dna: true,
            allow_ambiguity: false,
            ignore_softmask: false,
            sequence_starts: seq_data.start_positions,
            sequence_names: seq_data.sequence_names,
            num_partitions: 1,
            seed_mask: Some("101".to_string()),
            random_seed: 0,
        };

        // 7 $
        // 6 A$
        // 5 AA$
        // 4 AAA$
        // 2 ACAAA$
        // 0 AAACAAA$
        // 1 AACAAA$
        // 3 CAAA$

        let res = SufrBuilder::<u32>::new(builder_args);
        assert!(res.is_ok());

        let builder = res.unwrap();
        assert_eq!(builder.num_suffixes, 8);

        let mut sufr_file: SufrFile<u32> = SufrFile::read(&outpath, false)?;
        let suffix_array: Vec<_> = sufr_file.suffix_array_file.iter().collect();
        assert_eq!(suffix_array.len(), 8);
        assert_eq!(suffix_array, vec![7, 6, 5, 4, 2, 0, 1, 3]);

        Ok(())
    }

    #[test]
    fn test_spaced_seeds_2() -> Result<()> {
        let seq_file = Path::new("../data/inputs/mostlya2.fa");
        let sequence_delimiter = b'N';
        let seq_data = read_sequence_file(seq_file, sequence_delimiter)?;
        let outfile = NamedTempFile::new()?;
        let outpath = outfile.path().to_string_lossy().to_string();
        let builder_args = SufrBuilderArgs {
            text: seq_data.seq,
            low_memory: true,
            path: Some(outpath.clone()),
            max_query_len: None,
            is_dna: true,
            allow_ambiguity: false,
            ignore_softmask: false,
            sequence_starts: seq_data.start_positions,
            sequence_names: seq_data.sequence_names,
            num_partitions: 1,
            seed_mask: Some("11011".to_string()),
            random_seed: 0,
        };

        //  0 16 $
        //  1 13 AAC$
        //  2  5 AACAAACAAAC$
        //  3  1 AACAAACAAACAAAC$
        //  4  9 AACAAAC$
        //  5 12 AAAC$
        //  6  8 AAACAAAC$
        //  7  4 AAACAAACAAAC$
        //  8  0 AAACAAACAAACAAAC$
        //  9 14 AC$
        // 10 10 ACAAAC$
        // 11  6 ACAAACAAAC$
        // 12  2 ACAAACAAACAAAC$
        // 13 15 C$
        // 14 11 CAAAC$
        // 15  7 CAAACAAAC$
        // 16  3 CAAACAAACAAAC$

        let res = SufrBuilder::<u32>::new(builder_args);
        assert!(res.is_ok());

        let builder = res.unwrap();
        assert_eq!(builder.num_suffixes, 17);

        let mut sufr_file = SufrFile::<u32>::read(&outpath, false)?;
        let suffix_array: Vec<_> = sufr_file.suffix_array_file.iter().collect();
        assert_eq!(suffix_array.len(), 17);
        assert_eq!(
            suffix_array,
            vec![16, 13, 9, 5, 1, 12, 8, 4, 0, 14, 10, 6, 2, 15, 11, 7, 3],
        );

        Ok(())
    }

    // --------------------------------------------------
    #[test]
    fn test_spaced_seeds_3() -> Result<()> {
        let seq_file = Path::new("../data/inputs/spaced_input.fa");
        let sequence_delimiter = b'N';
        let seq_data = read_sequence_file(seq_file, sequence_delimiter)?;
        let outfile = NamedTempFile::new()?;
        let outpath = outfile.path().to_string_lossy().to_string();
        let builder_args = SufrBuilderArgs {
            text: seq_data.seq,
            low_memory: false,
            path: Some(outpath.clone()),
            max_query_len: None,
            is_dna: true,
            allow_ambiguity: false,
            ignore_softmask: false,
            sequence_starts: seq_data.start_positions,
            sequence_names: seq_data.sequence_names,
            num_partitions: 1,
            seed_mask: Some("11000111".to_string()),
            random_seed: 0,
        };

        let res = SufrBuilder::<u32>::new(builder_args);
        assert!(res.is_ok());

        let builder = res.unwrap();
        assert_eq!(builder.num_suffixes, 43);

        let mut sufr_file: SufrFile<u32> = SufrFile::read(&outpath, false)?;
        let suffix_array: Vec<_> = sufr_file.suffix_array_file.iter().collect();
        assert_eq!(suffix_array.len(), 43);
        assert_eq!(
            suffix_array,
            vec![
                42, 18, 12, 0, 32, 29, 13, 23, 21, 6, 40, 1, 33, 19, 30, 10, 28, 9, 17,
                14, 4, 26, 39, 22, 25, 38, 24, 35, 7, 36, 15, 41, 5, 20, 31, 11, 27, 8,
                16, 3, 37, 34, 2
            ],
        );

        Ok(())
    }
}
