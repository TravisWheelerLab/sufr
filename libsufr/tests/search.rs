//! Search results on inputs small enough to read, every expectation
//! worked out from the text by hand.

use anyhow::Result;
use libsufr::{
    Alphabet, Context, Hits, OpenOptions, SeedMask, Storage, SuffixArray,
    SuffixArrayBuilder,
};
use std::fs;
use std::sync::OnceLock;
use tempfile::{NamedTempFile, TempDir};

/// The DNA alphabet with nothing allowed.
const DNA: Alphabet = Alphabet::Dna {
    allow_ambiguity: false,
    ignore_softmask: false,
};

/// Every way a file can be opened: text and suffixes each in memory or
/// on disk, with and without a cap at search time. A cap with the
/// suffixes in memory is the subsampled path.
const MODES: [(Storage, Storage, Option<usize>); 8] = [
    (Storage::Memory, Storage::Memory, None),
    (Storage::Memory, Storage::Disk, None),
    (Storage::Disk, Storage::Memory, None),
    (Storage::Disk, Storage::Disk, None),
    (Storage::Memory, Storage::Memory, Some(2)),
    (Storage::Memory, Storage::Disk, Some(2)),
    (Storage::Disk, Storage::Memory, Some(2)),
    (Storage::Disk, Storage::Disk, Some(2)),
];

/// Build `fasta` under data/inputs, or a FASTA string, to a temp file.
fn build(fasta: &str, alphabet: Alphabet, mask: Option<&str>) -> Result<NamedTempFile> {
    let input;
    let path = if fasta.starts_with('>') {
        input = NamedTempFile::new()?;
        fs::write(input.path(), fasta)?;
        input.path().to_path_buf()
    } else {
        format!("../data/inputs/{fasta}.fa").into()
    };
    let mut builder = SuffixArrayBuilder::from_sequence_file(path, b'%')?
        .with_alphabet(alphabet)
        .with_write_lcp(true);
    if let Some(mask) = mask {
        builder = builder.with_seed_mask(SeedMask::new(mask)?);
    }
    let outfile = NamedTempFile::new()?;
    builder.build_to_file(outfile.path(), false)?;
    Ok(outfile)
}

fn open(
    file: &NamedTempFile,
    text: Storage,
    suffixes: Storage,
    max_query_len: Option<usize>,
) -> Result<SuffixArray> {
    // TODO: a subsampled search writes a cache under ~/.sufr; HOME is
    //       a temp dir for this process until sufr-22b9c3 decides
    //       where, or whether, that cache lives
    static HOME: OnceLock<TempDir> = OnceLock::new();
    HOME.get_or_init(|| {
        let home = TempDir::new().expect("a temp dir for HOME");
        std::env::set_var("HOME", home.path());
        home
    });
    SuffixArray::open_with(
        &file.path().to_string_lossy(),
        &OpenOptions {
            text,
            suffixes,
            max_query_len,
        },
    )
}

/// Hits as (sequence name, position in it, position in the text), in
/// the order locate gives them.
fn located(hits: &Hits) -> Result<Vec<(String, usize, usize)>> {
    Ok(hits
        .locate()?
        .iter()
        .map(|l| (l.sequence_name.to_string(), l.position, l.suffix))
        .collect())
}

/// Extractions as (sequence name, window, offset of the hit), sorted,
/// since extract gives them in rank order.
fn extracted(
    hits: &Hits,
    context: Context,
) -> Result<Vec<(String, (usize, usize), usize)>> {
    let mut windows: Vec<_> = hits
        .extract(context)?
        .iter()
        .map(|e| {
            (
                e.location.sequence_name.to_string(),
                (e.sequence_range.start, e.sequence_range.end),
                e.suffix_offset,
            )
        })
        .collect();
    windows.sort();
    Ok(windows)
}

// the text of 2.fa is ACGTACGT%ACGTACGT$: ABC at 0..8, DEF at 9..17

#[test]
fn every_mode_finds_the_same_hits() -> Result<()> {
    let file = build("2", DNA, None)?;
    for (text, suffixes, max_query_len) in MODES {
        let label = format!("{text:?} {suffixes:?} {max_query_len:?}");
        let sa = open(&file, text, suffixes, max_query_len)?;
        let hits = sa.search(b"ACGT")?;
        assert_eq!(hits.count(), 4, "{label}");
        assert_eq!(
            located(&hits)?,
            [
                ("ABC".to_string(), 0, 0),
                ("ABC".to_string(), 4, 4),
                ("DEF".to_string(), 0, 9),
                ("DEF".to_string(), 4, 13),
            ],
            "{label}"
        );
        let context = Context {
            before: 1,
            after: Some(2),
        };
        assert_eq!(
            extracted(&hits, context)?,
            [
                ("ABC".to_string(), (0, 2), 0),
                ("ABC".to_string(), (3, 6), 1),
                ("DEF".to_string(), (0, 2), 0),
                ("DEF".to_string(), (3, 6), 1),
            ],
            "{label}"
        );
    }
    Ok(())
}

// the text of 1.fa is ACGTNNACGT$, one sequence named 1

#[test]
fn query_kinds() -> Result<()> {
    let plain = build("1", DNA, None)?;
    let ambiguous = build(
        "1",
        Alphabet::Dna {
            allow_ambiguity: true,
            ignore_softmask: false,
        },
        None,
    )?;
    for (text, suffixes, _) in &MODES[..4] {
        let label = format!("{text:?} {suffixes:?}");
        let sa = open(&plain, *text, *suffixes, None)?;
        // sorts before the sentinel, the first suffix
        assert_eq!(sa.search(b"!")?.count(), 0, "{label}");
        assert_eq!(sa.search(b"$")?.count(), 1, "{label}");
        assert_eq!(sa.search(b"ACGTNNACGT$")?.suffixes()?, [0], "{label}");
        // matches the suffix T$ to the end of the text with a byte
        // left over
        assert_eq!(sa.search(b"ACGTNNACGT$A")?.count(), 0, "{label}");
        assert_eq!(sa.search(b"T$A")?.count(), 0, "{label}");
        // N-suffixes are not indexed without allow_ambiguity
        assert_eq!(sa.search(b"N")?.count(), 0, "{label}");
        // T at 3 and at the text's end, 9
        assert_eq!(sa.search(b"T")?.suffixes()?, [9, 3], "{label}");
        // the empty query is a prefix of every suffix
        assert_eq!(sa.search(b"")?.ranks(), 0..9, "{label}");

        let sa = open(&ambiguous, *text, *suffixes, None)?;
        // NNACGT$ at 4 and NACGT$ at 5
        assert_eq!(sa.search(b"N")?.suffixes()?, [5, 4], "{label}");
        assert_eq!(sa.search(b"")?.ranks(), 0..11, "{label}");
    }
    Ok(())
}

#[test]
fn a_whole_sequence_matches_at_its_start() -> Result<()> {
    let file = build("2", DNA, None)?;
    let sa = open(&file, Storage::Memory, Storage::Memory, None)?;
    assert_eq!(
        located(&sa.search(b"ACGTACGT")?)?,
        [("ABC".to_string(), 0, 0), ("DEF".to_string(), 0, 9)]
    );
    Ok(())
}

#[test]
fn the_cap_at_search_time_bounds_the_comparison() -> Result<()> {
    let file = build("2", DNA, None)?;
    // no cap: X matches nothing
    let sa = open(&file, Storage::Memory, Storage::Memory, None)?;
    assert_eq!(sa.search(b"ACGTX")?.count(), 0);
    // a cap shorter than the query: only AC is compared
    let sa = open(&file, Storage::Memory, Storage::Memory, Some(2))?;
    assert_eq!(sa.search(b"ACGTX")?.count(), 4);
    // a cap equal to the query
    let sa = open(&file, Storage::Memory, Storage::Memory, Some(4))?;
    assert_eq!(sa.search(b"ACGT")?.count(), 4);
    assert_eq!(sa.search(b"ACGTA")?.count(), 4);
    // a cap longer than the query changes nothing
    let sa = open(&file, Storage::Memory, Storage::Memory, Some(8))?;
    assert_eq!(sa.search(b"ACGT")?.count(), 4);
    assert_eq!(sa.search(b"ACGTX")?.count(), 0);
    Ok(())
}

#[test]
fn a_cap_of_zero_at_search_time_means_the_files_cap() -> Result<()> {
    // built on the first 2 bytes, so no search may compare more
    let outfile = NamedTempFile::new()?;
    SuffixArrayBuilder::from_sequence_file("../data/inputs/2.fa", b'%')?
        .with_alphabet(DNA)
        .with_max_query_len(2)
        .build_to_file(outfile.path(), false)?;
    for max_query_len in [None, Some(0), Some(9)] {
        let sa = open(&outfile, Storage::Memory, Storage::Memory, max_query_len)?;
        assert_eq!(sa.search(b"ACGTX")?.count(), 4, "{max_query_len:?}");
    }
    Ok(())
}

#[test]
fn a_seed_mask_skips_its_dont_care_positions() -> Result<()> {
    // mask 101 reads bytes 0 and 2 of the query and the suffix
    let file = build("2", DNA, Some("101"))?;
    for (text, suffixes, _) in &MODES[..4] {
        let label = format!("{text:?} {suffixes:?}");
        let sa = open(&file, *text, *suffixes, None)?;
        // A?G matches every ACGT
        assert_eq!(sa.search(b"AGGT")?.suffixes()?.len(), 4, "{label}");
        // A?C matches nothing
        assert_eq!(sa.search(b"AGCT")?.count(), 0, "{label}");
    }
    Ok(())
}

#[test]
fn extract_windows_stop_at_the_sequence() -> Result<()> {
    let file = build("2", DNA, None)?;
    let sa = open(&file, Storage::Memory, Storage::Memory, None)?;
    let hits = sa.search(b"ACGT")?;
    // the hit at DEF:4, where DEF is 8 bytes long
    let on_def = |context| -> Result<_> {
        Ok(extracted(&hits, context)?
            .into_iter()
            .filter(|(name, _, _)| name == "DEF")
            .map(|(_, window, offset)| (window, offset))
            .last())
    };
    let cases = [
        (100, None, ((0, 8), 4)),
        (0, Some(0), ((4, 4), 0)),
        (0, Some(10), ((4, 8), 0)),
        (0, None, ((4, 8), 0)),
        (1, Some(1), ((3, 5), 1)),
    ];
    for (before, after, expected) in cases {
        let context = Context { before, after };
        assert_eq!(on_def(context)?, Some(expected), "{context:?}");
    }
    Ok(())
}

#[test]
fn locate_orders_sequences_as_the_file_did() -> Result<()> {
    let file = build(">1\nACGT\n>2\nACGT\n>10\nACGT\n", DNA, None)?;
    let sa = open(&file, Storage::Memory, Storage::Memory, None)?;
    let names: Vec<String> = located(&sa.search(b"ACGT")?)?
        .into_iter()
        .map(|(name, _, _)| name)
        .collect();
    assert_eq!(names, ["1", "2", "10"]);
    Ok(())
}
