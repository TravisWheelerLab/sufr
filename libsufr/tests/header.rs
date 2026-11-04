use anyhow::Result;
use libsufr::{
    Alphabet, OpenOptions, Storage, SuffixArray, SuffixArrayBuilder, SuffixSortType,
};
use tempfile::NamedTempFile;

fn build_and_read(write_lcp: bool) -> Result<SuffixArray> {
    let outfile = NamedTempFile::new()?;
    let outpath = outfile.path().to_string_lossy().to_string();
    SuffixArrayBuilder::from_sequence_file("../data/inputs/1.fa", b'%')?
        .with_alphabet(Alphabet::Dna {
            allow_ambiguity: false,
            ignore_softmask: false,
        })
        .with_write_lcp(write_lcp)
        .build_to_file(&outpath)?;
    SuffixArray::open_with(
        &outpath,
        &OpenOptions {
            text: Storage::Memory,
            suffixes: Storage::Disk,
            max_query_len: None,
        },
    )
}
fn check_common(sufr: &SuffixArray) {
    let header = sufr.header();
    assert_eq!(header.version, 6);
    assert!(header.is_dna);
    assert!(!header.allow_ambiguity);
    assert!(!header.ignore_softmask);
    assert_eq!(header.text_len, 11);
    assert_eq!(header.num_suffixes, 9);
    assert_eq!(header.sequence_starts, [0]);
    assert_eq!(header.sequence_names, ["1"]);
    assert_eq!(header.sort_type, SuffixSortType::MaxQueryLen(0));
    assert_eq!(sufr.text(0..11).unwrap(), b"ACGTNNACGT$");
    assert_eq!(header.text_pos, 72);
    assert_eq!(header.suffix_array_pos, 83);
}

#[test]
fn header_with_lcp() -> Result<()> {
    let sufr = build_and_read(true)?;
    check_common(&sufr);
    assert!(sufr.lcps(0..1).is_ok());
    assert_eq!(sufr.header().lcp_pos, Some(119));
    Ok(())
}

#[test]
fn header_without_lcp() -> Result<()> {
    let sufr = build_and_read(false)?;
    check_common(&sufr);
    assert!(sufr.lcps(0..1).is_err());
    assert_eq!(sufr.header().lcp_pos, None);
    Ok(())
}
