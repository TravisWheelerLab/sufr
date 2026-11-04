//! Every small input built under every flag set, checked against a
//! reference that applies the rules to the FASTA bytes and never calls
//! the build.

use anyhow::Result;
use libsufr::{Alphabet, SeedMask, SuffixArray, SuffixArrayBuilder};
use std::borrow::Cow;
use std::fs;
use tempfile::NamedTempFile;

const SENTINEL: u8 = b'$';

/// One build of the table.
struct Case {
    input: &'static str,

    /// `allow_ambiguity` and `ignore_softmask` for a DNA build, or
    /// `None` for the byte alphabet.
    dna: Option<(bool, bool)>,

    delimiter: u8,

    /// Zero for no cap.
    max_query_len: usize,

    mask: Option<&'static str>,

    /// `None` for the builder's default.
    partitions: Option<usize>,
}

/// The text the build should index: the records joined by the
/// delimiter with the sentinel last, lowercase raised or turned to N.
fn reference_text(
    fasta: &[u8],
    delimiter: u8,
    ignore_softmask: bool,
) -> (Vec<u8>, Vec<usize>, Vec<String>) {
    let mut text = Vec::new();
    let mut starts = Vec::new();
    let mut names = Vec::new();
    for record in fasta.split(|&b| b == b'>').skip(1) {
        let mut lines = record.split(|&b| b == b'\n');
        let header = String::from_utf8_lossy(lines.next().unwrap_or(&[]));
        if !starts.is_empty() {
            text.push(delimiter);
        }
        starts.push(text.len());
        names.push(header.split_whitespace().next().unwrap_or("").to_string());
        for &b in lines.flatten() {
            if b == b'\r' {
                continue;
            }
            text.push(if !b.is_ascii_lowercase() {
                b
            } else if ignore_softmask {
                b'N'
            } else {
                b.to_ascii_uppercase()
            });
        }
    }
    text.push(SENTINEL);
    (text, starts, names)
}

/// The positions the build should index: without `allow_ambiguity` a
/// DNA build keeps only suffixes starting `ACGT` or the sentinel.
fn kept(text: &[u8], dna: Option<(bool, bool)>) -> Vec<usize> {
    let acgt_only = matches!(dna, Some((false, _)));
    (0..text.len())
        .filter(|&p| {
            !acgt_only || matches!(text[p], b'A' | b'C' | b'G' | b'T' | SENTINEL)
        })
        .collect()
}

/// What the suffix at `p` sorts by: its first `max_query_len` bytes, or
/// the bytes at the mask's positions that fall inside the text.
fn key<'a>(
    text: &'a [u8],
    p: usize,
    max_query_len: usize,
    mask: Option<&[usize]>,
) -> Cow<'a, [u8]> {
    match mask {
        Some(positions) => positions
            .iter()
            .map(|&offset| p + offset)
            .take_while(|&i| i < text.len())
            .map(|i| text[i])
            .collect(),
        None if max_query_len > 0 => {
            Cow::Borrowed(&text[p..(p + max_query_len).min(text.len())])
        }
        None => Cow::Borrowed(&text[p..]),
    }
}

/// The number of leading units two keys share.
fn shared(a: &[u8], b: &[u8]) -> usize {
    a.iter().zip(b).take_while(|(x, y)| x == y).count()
}

fn check(case: &Case) -> Result<()> {
    let Case {
        input,
        dna,
        delimiter,
        max_query_len,
        mask,
        partitions,
    } = *case;
    let label = format!(
        "{input} dna={dna:?} delimiter={} max_query_len={max_query_len} mask={mask:?} partitions={partitions:?}",
        delimiter as char
    );
    let path = format!("../data/inputs/{input}.fa");
    let mask = mask.map(SeedMask::new).transpose()?;
    let positions = mask.as_ref().map(|m| m.positions.clone());

    let ignore_softmask = matches!(dna, Some((_, true)));
    let (text, starts, names) =
        reference_text(&fs::read(&path)?, delimiter, ignore_softmask);
    let mut order = kept(&text, dna);
    order.sort_by(|&a, &b| {
        key(&text, a, max_query_len, positions.as_deref()).cmp(&key(
            &text,
            b,
            max_query_len,
            positions.as_deref(),
        ))
    });

    let alphabet = match dna {
        Some((allow_ambiguity, ignore_softmask)) => Alphabet::Dna {
            allow_ambiguity,
            ignore_softmask,
        },
        None => Alphabet::Bytes,
    };
    let mut builder = SuffixArrayBuilder::from_sequence_file(&path, delimiter)?
        .with_alphabet(alphabet)
        .with_write_lcp(true);
    if max_query_len > 0 {
        builder = builder.with_max_query_len(max_query_len);
    }
    if let Some(mask) = mask {
        builder = builder.with_seed_mask(mask);
    }
    if let Some(partitions) = partitions {
        builder = builder.with_num_partitions(partitions);
    }
    let outfile = NamedTempFile::new()?;
    let outpath = outfile.path().to_string_lossy().to_string();
    builder.build_to_file(&outpath, false)?;

    let built = SuffixArray::open(&outpath)?;
    let header = built.header();
    assert_eq!(header.text_len, text.len(), "{label}: text_len");
    assert_eq!(built.text(0..text.len())?, text, "{label}: text");
    assert_eq!(header.sequence_starts, starts, "{label}: sequence_starts");
    assert_eq!(header.sequence_names, names, "{label}: sequence_names");
    assert_eq!(header.num_suffixes, order.len(), "{label}: num_suffixes");

    let n = order.len();
    let suffixes = built.suffixes(0..n)?;
    let lcps = built.lcps(0..n)?;
    let mut group = 0;
    for rank in 0..n {
        let expected = key(&text, order[rank], max_query_len, positions.as_deref());
        let actual = key(&text, suffixes[rank], max_query_len, positions.as_deref());
        assert_eq!(actual, expected, "{label}: key at rank {rank}");

        let lcp = if rank == 0 {
            0
        } else {
            shared(
                &key(&text, order[rank - 1], max_query_len, positions.as_deref()),
                &expected,
            )
        };
        assert_eq!(lcps[rank], lcp, "{label}: lcp at rank {rank}");

        // the order inside a run of equal keys is the merge's, so
        // a run is compared as a set of positions
        let last_of_group = rank + 1 == n
            || key(&text, order[rank + 1], max_query_len, positions.as_deref())
                != expected;
        if last_of_group {
            let mut expected: Vec<usize> = order[group..=rank].to_vec();
            let mut actual: Vec<usize> = suffixes[group..=rank].to_vec();
            expected.sort_unstable();
            actual.sort_unstable();
            assert_eq!(
                actual, expected,
                "{label}: positions at ranks {group}..={rank}"
            );
            group = rank + 1;
        }
    }
    Ok(())
}

const DNA_INPUTS: [&str; 7] = [
    "1",
    "2",
    "3",
    "smol",
    "mostlya1",
    "mostlya2",
    "spaced_input",
];

fn dna(input: &'static str, allow_ambiguity: bool, ignore_softmask: bool) -> Case {
    Case {
        input,
        dna: Some((allow_ambiguity, ignore_softmask)),
        delimiter: b'%',
        max_query_len: 0,
        mask: None,
        partitions: None,
    }
}

fn bytes(input: &'static str) -> Case {
    Case {
        input,
        dna: None,
        delimiter: b'%',
        max_query_len: 0,
        mask: None,
        partitions: None,
    }
}

fn table() -> Vec<Case> {
    let mut cases = Vec::new();
    for input in DNA_INPUTS {
        for (allow_ambiguity, ignore_softmask) in
            [(false, false), (true, false), (false, true), (true, true)]
        {
            cases.push(dna(input, allow_ambiguity, ignore_softmask));
        }
        for max_query_len in [1, 3, 5] {
            for allow_ambiguity in [false, true] {
                cases.push(Case {
                    max_query_len,
                    ..dna(input, allow_ambiguity, false)
                });
            }
        }
    }
    for input in ["2", "3"] {
        cases.push(Case {
            delimiter: b'N',
            ..dna(input, false, false)
        });
    }
    for input in ["2", "smol", "mostlya1", "mostlya2", "spaced_input"] {
        for mask in ["101", "1101", "11011", "10111011"] {
            for allow_ambiguity in [false, true] {
                cases.push(Case {
                    mask: Some(mask),
                    ..dna(input, allow_ambiguity, false)
                });
            }
        }
    }
    for partitions in [1, 2, 3] {
        cases.push(Case {
            partitions: Some(partitions),
            ..dna("2", false, false)
        });
        cases.push(Case {
            partitions: Some(partitions),
            max_query_len: 3,
            ..dna("2", false, false)
        });
        cases.push(Case {
            partitions: Some(partitions),
            mask: Some("101"),
            ..dna("2", false, false)
        });
    }
    cases.push(bytes("abba"));
    cases.push(Case {
        max_query_len: 2,
        ..bytes("abba")
    });
    cases.push(bytes("uniprot"));
    cases.push(Case {
        mask: Some("10111011"),
        ..bytes("uniprot")
    });
    cases
}

#[test]
fn every_build_matches_the_reference() -> Result<()> {
    for case in table() {
        check(&case)?;
    }
    Ok(())
}
