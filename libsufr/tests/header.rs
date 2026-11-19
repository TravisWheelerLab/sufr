use anyhow::Result;
use libsufr::{OpenOptions, SoftMask, SourceKind, SuffixArray, SuffixArrayBuilder};
use tempfile::NamedTempFile;

fn build_and_read(write_lcp: bool) -> Result<SuffixArray> {
    let outfile = NamedTempFile::new()?;
    let outpath = outfile.path().to_string_lossy().to_string();
    SuffixArrayBuilder::from_sequence_file("../data/inputs/1.fa", b'%')?
        .with_starts(b"ACGT")
        .with_write_lcp(write_lcp)
        .build_to_file(&outpath)?;
    Ok(SuffixArray::open_with(
        &outpath,
        &OpenOptions {
            text: SourceKind::Memory,
            suffixes: SourceKind::Disk,
            max_query_len: None,
        },
    )?)
}
fn check_common(sufr: &SuffixArray) {
    let header = sufr.header();
    assert_eq!(header.version, 7);
    assert_eq!(header.starts.as_deref(), Some(b"ACGT".as_slice()));
    assert_eq!(header.soft_mask, SoftMask::Fold);
    assert_eq!(header.text_len, 11);
    assert_eq!(header.num_suffixes, 9);
    assert_eq!(header.sequence_starts, [0]);
    assert_eq!(header.sequence_names, ["1"]);
    assert_eq!(header.seed_mask(), None);
    assert_eq!(sufr.text(0..11).unwrap(), b"ACGTNNACGT$");
}

#[test]
fn header_with_lcp() -> Result<()> {
    let sufr = build_and_read(true)?;
    check_common(&sufr);
    assert!(sufr.lcps(0..1).is_ok());
    assert!(sufr.header().has_lcp());
    Ok(())
}

#[test]
fn header_without_lcp() -> Result<()> {
    let sufr = build_and_read(false)?;
    check_common(&sufr);
    assert!(sufr.lcps(0..1).is_err());
    assert!(!sufr.header().has_lcp());
    Ok(())
}
