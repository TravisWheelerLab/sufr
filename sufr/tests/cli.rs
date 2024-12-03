use anyhow::Result;
use assert_cmd::Command;
use libsufr::types::OUTFILE_VERSION;
use pretty_assertions::assert_eq;
use regex::Regex;
use std::fs;
use tempfile::NamedTempFile;

const PRG: &str = "sufr";
const SEQ1: &str = "../data/inputs/1.fa";
const SEQ2: &str = "../data/inputs/2.fa";
const SEQ3: &str = "../data/inputs/3.fa";
const SUFR1: &str = "../data/expected/1.sufr";
const SUFR2: &str = "../data/expected/2.sufr";

struct CountOptions {
    queries: Vec<String>,
    low_memory: bool,
}

struct CreateOptions {
    allow_ambiguity: bool,
    ignore_softmask: bool,
    sequence_delimiter: Option<char>,
}

struct ExtractOptions {
    positions: Vec<String>,
    prefix_len: Option<usize>,
    suffix_len: Option<usize>,
}

struct ListOptions {
    show_rank: bool,
    show_suffix: bool,
    show_lcp: bool,
    suffix_len: Option<usize>,
}

struct LocateOptions {
    queries: Vec<String>,
    absolute: bool,
}

// --------------------------------------------------
fn create(input_file: &str, expected_file: &str, opts: CreateOptions) -> Result<()> {
    let outfile = NamedTempFile::new()?;
    let outpath = &outfile.path().to_string_lossy();
    let mut args = vec![
        "create".to_string(),
        "--dna".to_string(),
        "-o".to_string(),
        outpath.to_string(),
        input_file.to_string(),
    ];

    if opts.allow_ambiguity {
        args.push("--allow-ambiguity".to_string());
    }

    if opts.ignore_softmask {
        args.push("--ignore-softmask".to_string());
    }

    if let Some(delim) = opts.sequence_delimiter {
        let mut tmp = vec!["--sequence-delimiter".to_string(), delim.to_string()];
        args.append(&mut tmp);
    }

    let output = Command::cargo_bin(PRG)?.args(args).output().expect("fail");

    assert!(output.status.success());
    assert!(outfile.path().exists());

    let actual = fs::read(outfile.path())?;
    let expected = fs::read(expected_file)?;

    assert_eq!(actual, expected);

    Ok(())
}

// --------------------------------------------------
fn check(filename: &str) -> Result<()> {
    let output = Command::cargo_bin(PRG)?
        .args(["check", filename])
        .output()
        .expect("fail");

    assert!(output.status.success());

    let stdout = String::from_utf8(output.stdout).expect("invalid UTF-8");

    assert!(stdout.contains("found 0 errors"));

    Ok(())
}

// --------------------------------------------------
#[test]
fn create_empty_dies() -> Result<()> {
    Command::cargo_bin(PRG)?
        .args(["create", "../data/expected/empty.sa"])
        .assert()
        .failure();
    Ok(())
}

// --------------------------------------------------
#[test]
fn create_seq1() -> Result<()> {
    create(
        SEQ1,
        "../data/expected/1.sufr",
        CreateOptions {
            allow_ambiguity: false,
            ignore_softmask: false,
            sequence_delimiter: None,
        },
    )
}

// --------------------------------------------------
#[test]
fn create_seq1_allow_ambiguity() -> Result<()> {
    create(
        SEQ1,
        "../data/expected/1n.sufr",
        CreateOptions {
            allow_ambiguity: true,
            ignore_softmask: false,
            sequence_delimiter: None,
        },
    )
}

// --------------------------------------------------
#[test]
fn create_seq2() -> Result<()> {
    create(
        SEQ2,
        "../data/expected/2.sufr",
        CreateOptions {
            allow_ambiguity: false,
            ignore_softmask: false,
            sequence_delimiter: None,
        },
    )
}

// --------------------------------------------------
#[test]
fn create_seq2_sequence_delimiter() -> Result<()> {
    create(
        SEQ2,
        "../data/expected/2d.sufr",
        CreateOptions {
            allow_ambiguity: false,
            ignore_softmask: false,
            sequence_delimiter: Some('N'),
        },
    )
}

// --------------------------------------------------
#[test]
fn create_seq2_allow_ambiguity() -> Result<()> {
    create(
        SEQ2,
        "../data/expected/2n.sufr",
        CreateOptions {
            allow_ambiguity: true,
            ignore_softmask: false,
            sequence_delimiter: None,
        },
    )
}

// --------------------------------------------------
#[test]
fn create_seq2_ignore_softmask() -> Result<()> {
    create(
        SEQ2,
        "../data/expected/2s.sufr",
        CreateOptions {
            allow_ambiguity: false,
            ignore_softmask: true,
            sequence_delimiter: None,
        },
    )
}

// --------------------------------------------------
#[test]
fn create_seq2_allow_ambiguity_ignore_softmask() -> Result<()> {
    create(
        SEQ2,
        "../data/expected/2ns.sufr",
        CreateOptions {
            allow_ambiguity: true,
            ignore_softmask: true,
            sequence_delimiter: None,
        },
    )
}

// --------------------------------------------------
#[test]
fn create_seq3() -> Result<()> {
    create(
        SEQ3,
        "../data/expected/3.sufr",
        CreateOptions {
            allow_ambiguity: false,
            ignore_softmask: false,
            sequence_delimiter: None,
        },
    )
}

// --------------------------------------------------
#[test]
fn check_seq1() -> Result<()> {
    check("../data/expected/1.sufr")
}

// --------------------------------------------------
fn count(
    filename: &str,
    opts: CountOptions,
    expected_stdout: &str,
    expected_stderr: Option<&str>,
) -> Result<()> {
    let mut args = vec!["count".to_string(), filename.to_string()];

    if opts.low_memory {
        args.push("-l".to_string());
    }

    for query in opts.queries {
        args.push(query.to_string());
    }

    let output = Command::cargo_bin(PRG)?.args(&args).output().expect("fail");
    assert!(output.status.success());

    let stdout = String::from_utf8(output.stdout).expect("invalid UTF-8");
    assert_eq!(stdout, expected_stdout);

    if let Some(expected) = expected_stderr {
        let stderr = String::from_utf8(output.stderr).expect("invalid UTF-8");
        assert_eq!(stderr, expected);
    }

    Ok(())
}

// --------------------------------------------------
#[test]
fn count_seq1() -> Result<()> {
    count(
        SUFR1,
        CountOptions {
            queries: vec!["AC".to_string(), "X".to_string(), "GT".to_string()],
            low_memory: false,
        },
        "AC 2\nGT 2\n",
        Some("X not found\n"),
    )
}

// --------------------------------------------------
fn extract(
    filename: &str,
    opts: ExtractOptions,
    expected: &str,
    error: Option<&str>,
) -> Result<()> {
    let mut args = vec!["extract".to_string(), filename.to_string()];

    if let Some(prefix_len) = opts.prefix_len {
        args.push("-p".to_string());
        args.push(prefix_len.to_string());
    }

    if let Some(suffix_len) = opts.suffix_len {
        args.push("-s".to_string());
        args.push(suffix_len.to_string());
    }

    for position in opts.positions {
        args.push(position.to_string());
    }
    dbg!(&args);

    let output = Command::cargo_bin(PRG)?.args(&args).output().expect("fail");
    assert!(output.status.success());

    let stdout = String::from_utf8(output.stdout).expect("invalid UTF-8");
    assert_eq!(stdout, expected);

    if let Some(err) = error {
        let stderr = String::from_utf8(output.stderr).expect("invalid UTF-8");
        assert_eq!(stderr, err);
    }

    Ok(())
}

// --------------------------------------------------
#[test]
fn extract_seq1() -> Result<()> {
    extract(
        SUFR1,
        ExtractOptions {
            positions: vec!["AC".to_string(), "GT".to_string(), "XX".to_string()],
            prefix_len: None,
            suffix_len: None,
        },
        &[
            ">1:6-11 AC 0",
            "ACGT$",
            ">1:0-11 AC 0",
            "ACGTNNACGT$",
            ">1:8-11 GT 0",
            "GT$",
            ">1:2-11 GT 0",
            "GTNNACGT$",
            "",
        ]
        .join("\n"),
        Some("XX not found\n"),
    )
}

// --------------------------------------------------
#[test]
fn extract_seq1_prefix_1() -> Result<()> {
    extract(
        SUFR1,
        ExtractOptions {
            positions: vec!["AC".to_string(), "GT".to_string()],
            prefix_len: Some(1),
            suffix_len: None,
        },
        &[
            ">1:5-11 AC 1",
            "NACGT$",
            ">1:0-11 AC 0",
            "ACGTNNACGT$",
            ">1:7-11 GT 1",
            "CGT$",
            ">1:1-11 GT 1",
            "CGTNNACGT$",
            "",
        ]
        .join("\n"),
        None,
    )
}

// --------------------------------------------------
#[test]
fn extract_seq1_suf_3() -> Result<()> {
    extract(
        SUFR1,
        ExtractOptions {
            positions: vec!["AC".to_string(), "GT".to_string()],
            prefix_len: None,
            suffix_len: Some(3),
        },
        &[
            ">1:6-9 AC 0",
            "ACG",
            ">1:0-3 AC 0",
            "ACG",
            ">1:8-11 GT 0",
            "GT$",
            ">1:2-5 GT 0",
            "GTN",
            "",
        ]
        .join("\n"),
        None,
    )
}

// --------------------------------------------------
#[test]
fn extract_seq1_pre_1_suf_3() -> Result<()> {
    extract(
        SUFR1,
        ExtractOptions {
            positions: vec!["AC".to_string(), "GT".to_string()],
            prefix_len: Some(1),
            suffix_len: Some(3),
        },
        &[
            ">1:5-9 AC 1",
            "NACG",
            ">1:0-3 AC 0",
            "ACG",
            ">1:7-11 GT 1",
            "CGT$",
            ">1:1-5 GT 1",
            "CGTN",
            "",
        ]
        .join("\n"),
        None,
    )
}

// --------------------------------------------------
fn list(filename: &str, opts: ListOptions, expected: &str) -> Result<()> {
    let mut args = vec!["list".to_string(), filename.to_string()];

    if opts.show_suffix {
        args.push("-s".to_string());
    }

    if opts.show_rank {
        args.push("-r".to_string());
    }

    if opts.show_lcp {
        args.push("-L".to_string());
    }

    if let Some(len) = opts.suffix_len {
        args.push("-l".to_string());
        args.push(len.to_string());
    }

    let output = Command::cargo_bin(PRG)?.args(&args).output().expect("fail");
    assert!(output.status.success());
    dbg!(&output);

    let stdout = String::from_utf8(output.stdout).expect("invalid UTF-8");
    assert_eq!(stdout, expected);

    Ok(())
}

// --------------------------------------------------
#[test]
fn list_sufr1_no_opts() -> Result<()> {
    list(
        SUFR1,
        ListOptions {
            show_lcp: false,
            show_suffix: false,
            show_rank: false,
            suffix_len: None,
        },
        &[
            "$",
            "ACGT$",
            "ACGTNNACGT$",
            "CGT$",
            "CGTNNACGT$",
            "GT$",
            "GTNNACGT$",
            "T$",
            "TNNACGT$",
            "",
        ]
        .join("\n"),
    )
}

// --------------------------------------------------
#[test]
fn list_sufr1_show_suffix() -> Result<()> {
    list(
        SUFR1,
        ListOptions {
            show_lcp: false,
            show_suffix: true,
            show_rank: false,
            suffix_len: None,
        },
        &[
            "10 $",
            " 6 ACGT$",
            " 0 ACGTNNACGT$",
            " 7 CGT$",
            " 1 CGTNNACGT$",
            " 8 GT$",
            " 2 GTNNACGT$",
            " 9 T$",
            " 3 TNNACGT$",
            "",
        ]
        .join("\n"),
    )
}

// --------------------------------------------------
#[test]
fn list_sufr1_show_lcp() -> Result<()> {
    list(
        SUFR1,
        ListOptions {
            show_lcp: true,
            show_suffix: false,
            show_rank: false,
            suffix_len: None,
        },
        &[
            " 0 $",
            " 0 ACGT$",
            " 4 ACGTNNACGT$",
            " 0 CGT$",
            " 3 CGTNNACGT$",
            " 0 GT$",
            " 2 GTNNACGT$",
            " 0 T$",
            " 1 TNNACGT$",
            "",
        ]
        .join("\n"),
    )
}

// --------------------------------------------------
#[test]
fn list_sufr1_show_rank_suffix_lcp() -> Result<()> {
    list(
        SUFR1,
        ListOptions {
            show_lcp: true,
            show_suffix: true,
            show_rank: true,
            suffix_len: None,
        },
        &[
            " 0 10  0 $",
            " 1  6  0 ACGT$",
            " 2  0  4 ACGTNNACGT$",
            " 3  7  0 CGT$",
            " 4  1  3 CGTNNACGT$",
            " 5  8  0 GT$",
            " 6  2  2 GTNNACGT$",
            " 7  9  0 T$",
            " 8  3  1 TNNACGT$",
            "",
        ]
        .join("\n"),
    )
}

// --------------------------------------------------
#[test]
fn list_sufr1_show_rank() -> Result<()> {
    list(
        SUFR1,
        ListOptions {
            show_lcp: false,
            show_suffix: false,
            show_rank: true,
            suffix_len: None,
        },
        &[
            " 0 $",
            " 1 ACGT$",
            " 2 ACGTNNACGT$",
            " 3 CGT$",
            " 4 CGTNNACGT$",
            " 5 GT$",
            " 6 GTNNACGT$",
            " 7 T$",
            " 8 TNNACGT$",
            "",
        ]
        .join("\n"),
    )
}
// --------------------------------------------------
#[test]
fn list_sufr1_len_3() -> Result<()> {
    list(
        SUFR1,
        ListOptions {
            show_lcp: false,
            show_suffix: false,
            show_rank: false,
            suffix_len: Some(3),
        },
        &[
            "$", "ACG", "ACG", "CGT", "CGT", "GT$", "GTN", "T$", "TNN", "",
        ]
        .join("\n"),
    )
}
// --------------------------------------------------
fn locate(filename: &str, opts: LocateOptions, expected: &str) -> Result<()> {
    let mut args = vec!["locate".to_string(), filename.to_string()];

    if opts.absolute {
        args.push("-a".to_string());
    }

    for query in opts.queries {
        args.push(query.to_string());
    }

    let output = Command::cargo_bin(PRG)?.args(&args).output().expect("fail");
    assert!(output.status.success());

    let stdout = String::from_utf8(output.stdout).expect("invalid UTF-8");
    assert_eq!(stdout, expected);

    Ok(())
}

// --------------------------------------------------
#[test]
fn locate_seq1_relative() -> Result<()> {
    locate(
        SUFR2,
        LocateOptions {
            queries: vec!["AC".to_string(), "GT".to_string()],
            absolute: false,
        },
        "AC\nABC 0,4\nDEF 0,4\n//\nGT\nABC 2,6\nDEF 2,6\n//\n",
    )
}

// --------------------------------------------------
#[test]
fn locate_seq1_absolute() -> Result<()> {
    locate(
        SUFR2,
        LocateOptions {
            queries: vec!["AC".to_string(), "GT".to_string()],
            absolute: true,
        },
        "AC 13 4 9 0\nGT 15 6 11 2\n",
    )
}

// --------------------------------------------------
fn summarize(filename: &str, expected: Vec<(&str, &str)>) -> Result<()> {
    let args = vec!["summarize".to_string(), filename.to_string()];

    let output = Command::cargo_bin(PRG)?.args(&args).output().expect("fail");
    assert!(output.status.success());

    let stdout = String::from_utf8(output.stdout).expect("invalid UTF-8");

    for (row_name, value) in expected {
        let pattern = format!(r#"[|] {}\s+[|] ([^|]+)"#, row_name);
        let regex = Regex::new(&pattern).unwrap();
        if let Some(caps) = regex.captures(&stdout) {
            assert_eq!(caps.get(1).unwrap().as_str().trim(), value);
        }
    }

    Ok(())
}

// --------------------------------------------------
#[test]
fn summarize_sufr1() -> Result<()> {
    summarize(
        SUFR1,
        vec![
            ("File Size", "172 bytes"),
            ("File Version", &OUTFILE_VERSION.to_string()),
            ("DNA", "true"),
            ("Allow Ambiguity", "false"),
            ("Ignore Softmask", "false"),
            ("Text Length", "11"),
            ("Num Suffixes", "9"),
            ("Max query len", "0"),
            ("Seed mask", "None"),
            ("Num sequences", "1"),
            ("Sequence starts", "0"),
            ("Headers", "1"),
        ],
    )
}
