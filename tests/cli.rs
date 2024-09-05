use anyhow::Result;
use assert_cmd::Command;
use std::fs;
use tempfile::NamedTempFile;

const PRG: &str = "sufr";
const EMPTY: &str = "tests/inputs/empty";
const SEQ1: &str = "tests/inputs/seq1.txt";
const SEQ1_SA: &str = "tests/inputs/seq1.sa";

// --------------------------------------------------
fn create(args: &[&str], expected_file: &str) -> Result<()> {
    let outfile = NamedTempFile::new()?;
    let outpath = &outfile.path().to_str().unwrap();
    let mut run_args = vec!["create", "-o", outpath];
    run_args.extend_from_slice(args);
    let output = Command::cargo_bin(PRG)?
        .args(&run_args)
        .output()
        .expect("fail");

    assert!(output.status.success());

    assert!(outfile.path().exists());

    let actual = fs::read(outfile.path())?;
    let expected = fs::read(expected_file)?;

    assert_eq!(actual, expected);

    Ok(())
}

// --------------------------------------------------
fn check(seq_file: &str, sa_file: &str) -> Result<()> {
    let args = vec!["check", "-s", seq_file, "-a", sa_file];
    let output = Command::cargo_bin(PRG)?.args(&args).output().expect("fail");

    assert!(output.status.success());

    let stdout = String::from_utf8(output.stdout).expect("invalid UTF-8");

    assert!(stdout.starts_with("Found 0 errors"));

    Ok(())
}

// --------------------------------------------------
fn read(
    seq_file: &str,
    sa_file: &str,
    extract: &str,
    expected_file: &str,
) -> Result<()> {
    let outfile = NamedTempFile::new()?;
    let outpath = &outfile.path().to_str().unwrap();
    let args = vec![
        "read", "-s", seq_file, "-a", sa_file, "-e", extract, "-o", outpath,
    ];
    let output = Command::cargo_bin(PRG)?.args(&args).output().expect("fail");

    assert!(output.status.success());
    assert!(outfile.path().exists());

    let actual = fs::read_to_string(outfile.path())?;
    let expected = fs::read_to_string(expected_file)?;

    assert_eq!(actual, expected);

    Ok(())
}

// --------------------------------------------------
#[test]
fn create_empty() -> Result<()> {
    create(&[EMPTY], "tests/expected/empty.sa")
}

// --------------------------------------------------
#[test]
fn create_seq1() -> Result<()> {
    create(&[SEQ1], "tests/expected/seq1.sa")
}

// --------------------------------------------------
#[test]
fn read_input1_1_10() -> Result<()> {
    read(SEQ1, SEQ1_SA, "1-10", "tests/expected/seq1.read.1-10.out")
}

// --------------------------------------------------
#[test]
fn read_input1_1_5() -> Result<()> {
    read(SEQ1, SEQ1_SA, "1-5", "tests/expected/seq1.read.1-5.out")
}

// --------------------------------------------------
#[test]
fn read_input1_5_10() -> Result<()> {
    read(SEQ1, SEQ1_SA, "5-10", "tests/expected/seq1.read.5-10.out")
}

// --------------------------------------------------
#[test]
fn check_input1() -> Result<()> {
    check(SEQ1, SEQ1_SA)
}
