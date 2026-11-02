use anyhow::Result;
use libsufr::{sufr_file::SufrFile, types::SuffixSortType, SuffixArrayBuilder};
use tempfile::NamedTempFile;

// --------------------------------------------------
fn build_and_read(write_lcp: bool) -> Result<SufrFile<u32>> {
    let outfile = NamedTempFile::new()?;
    let outpath = outfile.path().to_string_lossy().to_string();
    SuffixArrayBuilder::from_sequence_file("../data/inputs/1.fa", b'%')?
        .with_dna_options(true, false, false)
        .with_write_lcp(write_lcp)
        .build_to_file(&outpath, true)?;
    SufrFile::read(&outpath, false)
}

// --------------------------------------------------
fn check_common(sufr: &SufrFile<u32>) {
    assert_eq!(sufr.version, 6);
    assert!(sufr.is_dna);
    assert!(!sufr.allow_ambiguity);
    assert!(!sufr.ignore_softmask);
    assert_eq!(sufr.text_len, 11);
    assert_eq!(sufr.len_suffixes, 9);
    assert_eq!(sufr.num_sequences, 1);
    assert_eq!(sufr.sequence_starts, [0]);
    assert_eq!(sufr.sequence_names, ["1"]);
    assert_eq!(sufr.sort_type, SuffixSortType::MaxQueryLen(0));
    assert_eq!(sufr.text, b"ACGTNNACGT$");
    assert_eq!(sufr.text_pos, 72);
    assert_eq!(sufr.suffix_array_pos, 83);
}

// --------------------------------------------------
#[test]
fn header_with_lcp() -> Result<()> {
    let sufr = build_and_read(true)?;
    check_common(&sufr);
    assert!(sufr.has_lcp);
    assert_eq!(sufr.lcp_pos, 119);
    Ok(())
}

// --------------------------------------------------
#[test]
fn header_without_lcp() -> Result<()> {
    let sufr = build_and_read(false)?;
    check_common(&sufr);
    assert!(!sufr.has_lcp);
    assert_eq!(sufr.lcp_pos, 0);
    Ok(())
}
