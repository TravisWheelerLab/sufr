use anyhow::Result;
use assert_cmd::Command;
use pretty_assertions::assert_eq;
use std::fs;
use tempfile::NamedTempFile;

const PRG: &str = "sufr";
const SEQ1: &str = "tests/inputs/1.fa";
const SEQ2: &str = "tests/inputs/2.fa";
const SEQ3: &str = "tests/inputs/3.fa";
const SUFR2: &str = "tests/expected/2.sufr";

struct CreateOptions {
    allow_ambiguity: bool,
    ignore_softmask: bool,
    sequence_delimiter: Option<char>,
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
        .args(["create", "tests/expected/empty.sa"])
        .assert()
        .failure();
    Ok(())
}

// --------------------------------------------------
#[test]
fn create_seq1() -> Result<()> {
    create(
        SEQ1,
        "tests/expected/1.sufr",
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
        "tests/expected/1n.sufr",
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
        "tests/expected/2.sufr",
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
        "tests/expected/2d.sufr",
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
        "tests/expected/2n.sufr",
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
        "tests/expected/2s.sufr",
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
        "tests/expected/2ns.sufr",
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
        "tests/expected/3.sufr",
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
    check("tests/expected/1.sufr")
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
