use anyhow::Result;
use libsufr::{Builder, Index, Sort, Storage};
use tempfile::NamedTempFile;

fn build_and_read(write_lcp: bool) -> Result<Index<u32>> {
    let outfile = NamedTempFile::new()?;
    let outpath = outfile.path().to_string_lossy().to_string();
    Builder::from_sequence_file("../data/inputs/1.fa", b'%')?
        .dna(true, false, false)
        .write_lcp(write_lcp)
        .write(&outpath)?;
    Index::open(&outpath, Storage::Memory, Storage::Disk, None)
}

fn check_common(sufr: &Index<u32>) {
    let header = &sufr.header;
    assert_eq!(header.version, 6);
    assert!(header.is_dna);
    assert!(!header.allow_ambiguity);
    assert!(!header.ignore_softmask);
    assert_eq!(header.text_len, 11);
    assert_eq!(header.num_suffixes, 9);
    assert_eq!(header.sequence_starts, [0]);
    assert_eq!(header.sequence_names, ["1"]);
    assert_eq!(header.sort, Sort::MaxQueryLen(0));
    assert_eq!(sufr.text_range(0..11).unwrap(), b"ACGTNNACGT$");
    assert_eq!(header.text_pos, 72);
    assert_eq!(header.suffix_array_pos, 83);
}

#[test]
fn header_with_lcp() -> Result<()> {
    let sufr = build_and_read(true)?;
    check_common(&sufr);
    assert!(sufr.lcp_file.is_some());
    assert_eq!(sufr.header.lcp_pos, Some(119));
    Ok(())
}

#[test]
fn header_without_lcp() -> Result<()> {
    let sufr = build_and_read(false)?;
    check_common(&sufr);
    assert!(sufr.lcp_file.is_none());
    assert_eq!(sufr.header.lcp_pos, None);
    Ok(())
}
