pub mod sufr_builder;
pub mod sufr_file;
pub mod sufr_search;
pub mod file_access;
pub mod types;
pub mod util;

// --------------------------------------------------
#[cfg(test)]
mod tests {
    use super::{
        sufr_builder::{SufrBuilder, SufrBuilderArgs},
        sufr_file::SufrFile,
        types::OUTFILE_VERSION,
        util::read_sequence_file,
    };
    use anyhow::Result;
    use tempfile::NamedTempFile;

    #[test]
    fn test_write_read_suffix_file_32() -> Result<()> {
        let seq_file = "tests/inputs/2.fa";
        let sequence_delimiter = b'N';
        let seq_data = read_sequence_file(seq_file, sequence_delimiter)?;
        let args = SufrBuilderArgs {
            text: seq_data.seq,
            max_query_len: None,
            is_dna: true,
            allow_ambiguity: false,
            ignore_softmask: false,
            sequence_starts: seq_data.start_positions.into_iter().collect(),
            headers: seq_data.headers,
            num_partitions: 2,
            sequence_delimiter,
            seed_mask: None,
        };
        let sufr_builder: SufrBuilder<u32> = SufrBuilder::new(args)?;
        let outfile = NamedTempFile::new()?;
        let outpath = &outfile.path().to_str().unwrap();
        let res = sufr_builder.write(outpath);
        assert!(res.is_ok());
        assert!(outfile.path().exists());

        let res: Result<SufrFile<u32>> = SufrFile::read(outpath);
        assert!(res.is_ok());

        let mut sufr_file = res.unwrap();
        assert_eq!(sufr_file.version, OUTFILE_VERSION);
        assert!(sufr_file.is_dna);
        assert_eq!(sufr_file.text_len, 18);
        assert_eq!(sufr_file.num_sequences, 2);
        assert_eq!(sufr_file.sequence_starts, [0, 9]);
        assert_eq!(sufr_file.headers, ["ABC", "DEF"]);
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
        let seq_file = "tests/inputs/1.fa";
        let sequence_delimiter = b'N';
        let seq_data = read_sequence_file(seq_file, sequence_delimiter)?;
        let args = SufrBuilderArgs {
            text: seq_data.seq,
            max_query_len: None,
            is_dna: true,
            allow_ambiguity: true,
            ignore_softmask: false,
            sequence_starts: seq_data.start_positions.into_iter().collect(),
            headers: seq_data.headers,
            num_partitions: 2,
            sequence_delimiter,
            seed_mask: None,
        };

        let suffix_array: SufrBuilder<u64> = SufrBuilder::new(args)?;
        let sorted_sa = [10, 6, 0, 7, 1, 8, 2, 5, 4, 9, 3];
        let lcp = [0, 0, 4, 0, 3, 0, 2, 0, 1, 0, 1];
        let outfile = NamedTempFile::new()?;
        let outpath = &outfile.path().to_str().unwrap();
        let res = suffix_array.write(outpath);
        assert!(res.is_ok());
        assert!(outfile.path().exists());

        let res: Result<SufrFile<u64>> = SufrFile::read(outpath);
        assert!(res.is_ok());

        let mut sufr_file = res.unwrap();
        assert_eq!(sufr_file.version, OUTFILE_VERSION);
        assert!(sufr_file.is_dna);
        assert_eq!(sufr_file.text_len, 11);
        assert_eq!(sufr_file.num_sequences, 1);
        assert_eq!(sufr_file.sequence_starts, [0]);
        assert_eq!(sufr_file.headers, ["1"]);
        assert_eq!(sufr_file.text, b"ACGTNNACGT$");

        let file_sa: Vec<_> = sufr_file.suffix_array_file.iter().collect();
        assert_eq!(file_sa, sorted_sa);
        let file_lcp: Vec<_> = sufr_file.lcp_file.iter().collect();
        assert_eq!(file_lcp, lcp);
        Ok(())
    }

    #[test]
    fn test_subsample_suffix_array() -> Result<()> {
        let seq_file = "tests/inputs/smol.fa";
        let sequence_delimiter = b'N';
        let seq_data = read_sequence_file(seq_file, sequence_delimiter)?;
        let builder_args = SufrBuilderArgs {
            text: seq_data.seq,
            max_query_len: None,
            is_dna: true,
            allow_ambiguity: false,
            ignore_softmask: false,
            sequence_starts: seq_data.start_positions,
            headers: seq_data.headers,
            num_partitions: 2,
            sequence_delimiter,
            seed_mask: None,
        };
        let sufr: SufrBuilder<u32> = SufrBuilder::new(builder_args)?;
        let outfile = NamedTempFile::new()?;
        let outpath = &outfile.path().to_str().unwrap();
        let res = sufr.write(outpath);
        assert!(res.is_ok());
        assert_eq!(sufr.num_suffixes, 364);

        let mut sufr_file: SufrFile<u32> = SufrFile::read(outpath)?;
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
}
