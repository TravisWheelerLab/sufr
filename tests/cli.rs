use anyhow::Result;
use assert_cmd::Command;
use std::fs;
use tempfile::NamedTempFile;

const PRG: &str = "sufr";
const SEQ1: &str = "tests/inputs/1.fa";
const SEQ2: &str = "tests/inputs/2.fa";
const SEQ3: &str = "tests/inputs/3.fa";

// --------------------------------------------------
fn create(input_file: &str, expected_file: &str) -> Result<()> {
    let outfile = NamedTempFile::new()?;
    let outpath = &outfile.path().to_str().unwrap();
    let output = Command::cargo_bin(PRG)?
        .args(vec!["create", "--dna", "-o", outpath, input_file])
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
fn check(filename: &str) -> Result<()> {
    let output = Command::cargo_bin(PRG)?
        .args(&["check", filename])
        .output()
        .expect("fail");

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

    let actual = fs::read(outfile.path())?;
    let expected = fs::read(expected_file)?;

    assert_eq!(actual, expected);

    Ok(())
}

// --------------------------------------------------
#[test]
fn create_empty_dies() -> Result<()> {
    Command::cargo_bin(PRG)?
        .args(&["create", "tests/expected/empty.sa"])
        .assert()
        .failure();
    Ok(())
}

// --------------------------------------------------
#[test]
fn create_seq1() -> Result<()> {
    create(SEQ1, "tests/expected/1.sufr")
}

// --------------------------------------------------
#[test]
fn create_seq2() -> Result<()> {
    create(SEQ2, "tests/expected/2.sufr")
}

// --------------------------------------------------
#[test]
fn create_seq3() -> Result<()> {
    create(SEQ3, "tests/expected/3.sufr")
}

// --------------------------------------------------
#[test]
fn check_seq1() -> Result<()> {
    check("tests/expected/1.sufr")
}

// --------------------------------------------------
//#[test]
//fn read_input1_1_10() -> Result<()> {
//    read(SEQ1, SEQ1_SA, "1-10", "tests/expected/seq1.read.1-10.out")
//}

// --------------------------------------------------
//#[test]
//fn read_input1_1_5() -> Result<()> {
//    read(SEQ1, SEQ1_SA, "1-5", "tests/expected/seq1.read.1-5.out")
//}

// --------------------------------------------------
//#[test]
//fn read_input1_5_10() -> Result<()> {
//    read(SEQ1, SEQ1_SA, "5-10", "tests/expected/seq1.read.5-10.out")
//}

// --------------------------------------------------
//#[test]
//fn check_input1() -> Result<()> {
//    check(SEQ1, SEQ1_SA)
//}
