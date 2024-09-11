mod suffix_array;

use anyhow::{anyhow, bail, Result};
use clap::{builder::PossibleValue, Parser, ValueEnum};
use log::{debug, info};
use regex::Regex;
use std::{
    fs::File,
    io::{self, BufRead, BufReader, BufWriter, Read, Write},
    mem,
    num::NonZeroUsize,
    ops::Range,
    slice,
    time::Instant,
};
use substring::Substring;
use suffix_array::SuffixArray;

// --------------------------------------------------
#[derive(Parser, Debug)]
#[command(arg_required_else_help = true)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Option<Command>,
}

#[derive(Parser, Debug)]
pub enum Command {
    /// Create suffix array
    #[clap(alias = "ch")]
    Check(CheckArgs),

    /// Create suffix array
    #[clap(alias = "cr")]
    Create(CreateArgs),

    /// Read suffix array and extract sequences
    #[clap(alias = "re")]
    Read(ReadArgs),
}

#[derive(Debug, Parser)]
#[command(author, version, about)]
pub struct CheckArgs {
    /// Suffix array file
    #[arg(short, long, value_name = "SA")]
    pub array: String,

    /// Sequence file
    #[arg(short, long, value_name = "SEQ")]
    pub sequence: String,
}

#[derive(Debug, Parser)]
#[command(author, version, about)]
pub struct CreateArgs {
    /// Input file
    #[arg(value_name = "INPUT")]
    pub input: String,

    /// Subproblem count
    #[arg(short, long, value_name = "SUBPROBLEMS", default_value = "16")]
    pub subproblem_count: usize,

    /// Max context
    #[arg(short, long, value_name = "CONTEXT")]
    pub max_context: Option<usize>,

    /// Number of threads
    #[arg(short, long, value_name = "THREADS", default_value = "16")]
    pub threads: Option<usize>,

    /// Output file
    #[arg(short, long, value_name = "OUTPUT", default_value = "sufr.sa")]
    pub output: String,

    /// Ignore sequences starting with N
    #[arg(short, long)]
    pub ignore_start_n: bool,

    /// Verify order
    #[arg(short, long)]
    pub check: bool,

    /// Log level
    #[arg(short, long)]
    pub log: Option<LogLevel>,

    /// Log file
    #[arg(long)]
    pub log_file: Option<String>,
}

/// Read suffix array
#[derive(Debug, Parser)]
#[command(author, version, about)]
pub struct ReadArgs {
    /// Suffix array file
    #[arg(short, long, value_name = "SA")]
    pub array: String,

    /// Sequence file
    #[arg(short, long, value_name = "SEQ")]
    pub sequence: String,

    /// Maximum length of sequence
    #[arg(short, long, value_name = "MAX", default_value = "0")]
    pub max_len: usize,

    /// Extract positions
    #[arg(short, long, value_name = "EXTRACT", default_value = "1")]
    pub extract: String,

    /// Number output
    #[arg(short, long)]
    pub number: bool,

    /// Output
    #[arg(short, long)]
    pub output: Option<String>,
}

#[derive(Debug, Clone)]
pub enum LogLevel {
    Info,
    Debug,
}

impl ValueEnum for LogLevel {
    fn value_variants<'a>() -> &'a [Self] {
        &[LogLevel::Info, LogLevel::Debug]
    }

    fn to_possible_value<'a>(&self) -> Option<PossibleValue> {
        Some(match self {
            LogLevel::Info => PossibleValue::new("info"),
            LogLevel::Debug => PossibleValue::new("debug"),
        })
    }
}

type PositionList = Vec<Range<usize>>;

// --------------------------------------------------
pub fn check(args: &CheckArgs) -> Result<()> {
    let sa = read_suffix_array(&args.array)?;
    let sa_len = sa.len();
    let seq = SuffixArray::read_input(open(&args.sequence)?);

    if sa_len != seq.len() {
        bail!("SA len {sa_len} does not match sequence len {}", seq.len());
    }

    let start = Instant::now();
    let mut previous: Option<usize> = None;
    let mut num_errors = 0;
    for (i, &cur) in sa.iter().enumerate() {
        if let Some(p) = previous {
            if !is_less(&seq[p..sa_len], &seq[cur..sa_len]) {
                num_errors += 1;
                println!("POS {}", i + 1);
            }
        }
        previous = Some(cur);
    }

    println!(
        "Found {num_errors} error{} in {:?}.",
        if num_errors == 1 { "" } else { "s" },
        start.elapsed()
    );

    Ok(())
}

// --------------------------------------------------
fn is_less(s1: &[u8], s2: &[u8]) -> bool {
    let lcp = s1.iter().zip(s2).take_while(|(a, b)| a == b).count();
    match (s1.get(lcp), s2.get(lcp)) {
        (Some(a), Some(b)) => a < b,
        (None, Some(_)) => true,
        _ => false,
    }
}

#[test]
fn test_is_less() {
    assert!(is_less("$".as_bytes(), "A".as_bytes()));
    assert!(is_less("A".as_bytes(), "C".as_bytes()));
    assert!(!is_less("C".as_bytes(), "A".as_bytes()));
    assert!(is_less("A".as_bytes(), "AA".as_bytes()));
    assert!(!is_less("AA".as_bytes(), "A".as_bytes()));
    assert!(is_less("AA".as_bytes(), "AC".as_bytes()));
    assert!(!is_less("AC".as_bytes(), "AA".as_bytes()));
    assert!(is_less("AA".as_bytes(), "AAC".as_bytes()));
    assert!(!is_less("CA".as_bytes(), "AAC".as_bytes()));
}

// --------------------------------------------------
pub fn create(args: &CreateArgs) -> Result<()> {
    if let Some(num) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(num)
            .build_global()
            .unwrap();
    }

    env_logger::Builder::new()
        .filter_level(match args.log {
            Some(LogLevel::Debug) => log::LevelFilter::Debug,
            Some(LogLevel::Info) => log::LevelFilter::Info,
            _ => log::LevelFilter::Off,
        })
        .target(match args.log_file {
            // Optional log file, default to STDOUT
            Some(ref filename) => env_logger::Target::Pipe(Box::new(
                BufWriter::new(File::create(filename)?),
            )),
            _ => env_logger::Target::Stdout,
        })
        .init();

    let total_start = Instant::now();
    let start = Instant::now();
    let sa = SuffixArray::new(
        open(&args.input)?,
        args.max_context,
        args.ignore_start_n,
    );
    info!(
        "Read input of len {} with {} suffixes in {:?}",
        sa.len,
        sa.suffixes.len(),
        start.elapsed()
    );
    debug!("Raw input '{:?}'", sa.text);
    //dbg!(&sa);

    let start = Instant::now();
    let sufs = sa.sort_subarrays(args.subproblem_count);
    info!("Sorted subarrays in {:?}", start.elapsed());
    //dbg!(&sufs);
    debug!("sorted subarrays = {sufs:#?}");

    // Collect all the subarray SA/LCP structures
    let sub_suffixes: Vec<_> = sufs.iter().map(|t| t.0.clone()).collect();
    let sub_lcps: Vec<_> = sufs.iter().map(|t| t.1.clone()).collect();
    let sub_pivots: Vec<_> = sufs.iter().map(|t| t.2.clone()).collect();
    mem::drop(sufs);

    // Determine the pivot suffixes
    let start = Instant::now();
    let pivots = sa.select_pivots(sub_pivots, sub_suffixes.len());
    info!("Selected {} pivots in {:?}", pivots.len(), start.elapsed());
    debug!("pivots = {pivots:#?}");

    // Find the pivots suffixes in each of the sub_suffixes
    let start = Instant::now();
    let pivot_locs = sa.locate_pivots(&sub_suffixes, pivots);
    info!("Located pivot suffixes in {:?}", start.elapsed());
    debug!("pivot_locs = {pivot_locs:#?}");

    // Use the pivot locations to partition the SA/LCP subs
    let start = Instant::now();
    let (part_sas, part_lcps) =
        sa.partition_subarrays(&sub_suffixes, &sub_lcps, pivot_locs);
    info!("Partitioned subarrays in {:?}", start.elapsed());
    debug!("{part_sas:#?}");
    //mem::drop(sub_suffixes);
    //mem::drop(sub_lcps);
    //mem::drop(pivot_locs);

    // Merge the partitioned subs
    let start = Instant::now();
    let merged_sa = sa.merge_part_subs(part_sas, part_lcps);
    info!("Merged subarrays in {:?}", start.elapsed());
    debug!("{merged_sa:#?}");
    info!("Suffix generated in {:?}", total_start.elapsed());

    // Check
    if args.check {
        let start = Instant::now();
        let mut previous: Option<usize> = None;
        let mut num_errors = 0;
        for &cur in &merged_sa {
            if let Some(p) = previous {
                if !is_less(&sa.text[p..sa.len], &sa.text[cur..sa.len]) {
                    num_errors += 1;
                }
            }
            previous = Some(cur);
        }

        info!(
            "Checked order, found {num_errors} error{} in {:?}",
            if num_errors == 1 { "" } else { "s" },
            start.elapsed()
        );
    }

    let start = Instant::now();
    let mut out = File::create(&args.output)?;

    // Dynamically choose 32/64
    //let size = if sa.len < u32::MAX as usize {
    //    mem::size_of::<u32>()
    //} else {
    //    mem::size_of::<u64>()
    //};
    //let slice_u8: &[u8] = unsafe {
    //    slice::from_raw_parts(sa.as_ptr() as *const _, sa.len * size)
    //};

    // Write out suffix array length
    out.write(&usize_to_bytes(merged_sa.len()))?;

    // Write out suffix array as raw bytes
    let slice_u8: &[u8] = unsafe {
        slice::from_raw_parts(
            merged_sa.as_ptr() as *const _,
            merged_sa.len() * mem::size_of::<usize>(),
        )
    };
    out.write_all(slice_u8)?;
    info!("Wrote output file in {:?}", start.elapsed());

    println!("See output file '{}'", args.output);

    Ok(())
}

// --------------------------------------------------
fn usize_to_bytes(value: usize) -> Vec<u8> {
    // Determine the size of usize in bytes
    let size = std::mem::size_of::<usize>();

    // Create a vector to hold the bytes
    let mut bytes = Vec::with_capacity(size);

    // Convert usize to bytes
    for i in 0..size {
        bytes.push((value >> (i * 8)) as u8);
    }

    bytes
}

// --------------------------------------------------
fn open(filename: &str) -> Result<Box<dyn BufRead>> {
    match filename {
        "-" => Ok(Box::new(BufReader::new(io::stdin()))),
        _ => Ok(Box::new(BufReader::new(
            File::open(filename).map_err(|e| anyhow!("{filename}: {e}"))?,
        ))),
    }
}

// --------------------------------------------------
// Parse an index from a string representation of an integer.
// Ensures the number is non-zero.
// Ensures the number does not start with '+'.
// Returns an index, which is a non-negative integer that is
// one less than the number represented by the original input.
fn parse_index(input: &str) -> Result<usize> {
    let value_error = || anyhow!(r#"illegal list value: "{input}""#);
    input
        .starts_with('+')
        .then(|| Err(value_error()))
        .unwrap_or_else(|| {
            input
                .parse::<NonZeroUsize>()
                .map(|n| usize::from(n) - 1)
                .map_err(|_| value_error())
        })
}

// --------------------------------------------------
fn parse_pos(range: &str) -> Result<PositionList> {
    let range_re = Regex::new(r"^(\d+)-(\d+)$").unwrap();
    range
        .split(',')
        .map(|val| {
            parse_index(val).map(|n| n..n + 1).or_else(|e| {
                range_re.captures(val).ok_or(e).and_then(|captures| {
                    let n1 = parse_index(&captures[1])?;
                    let n2 = parse_index(&captures[2])?;
                    if n1 >= n2 {
                        bail!(
                            "First number in range ({}) \
                            must be lower than second number ({})",
                            n1 + 1,
                            n2 + 1
                        );
                    }
                    Ok(n1..n2 + 1)
                })
            })
        })
        .collect::<Result<_, _>>()
        .map_err(From::from)
}

// --------------------------------------------------
pub fn read(args: &ReadArgs) -> Result<()> {
    let start = Instant::now();
    let sa = read_suffix_array(&args.array)?;
    let sa_len = sa.len();
    info!("Read SA in {:?}", start.elapsed());

    let seq =
        String::from_utf8(SuffixArray::read_input(open(&args.sequence)?))?;

    if seq.len() != sa_len {
        bail!("SA len {sa_len} does not match sequence len {}", seq.len());
    }

    let positions: Vec<_> = parse_pos(&args.extract)?
        .into_iter()
        .flat_map(|r| r.collect::<Vec<_>>())
        .collect();

    let mut output: Box<dyn Write> = match &args.output {
        Some(out_name) => Box::new(File::create(out_name)?),
        _ => Box::new(io::stdout()),
    };

    for start in positions {
        if let Some(&pos) = sa.get(start) {
            //let pos = pos as usize;
            let end = if args.max_len > 0 {
                pos + args.max_len
            } else {
                sa_len
            };

            if args.number {
                writeln!(
                    output,
                    "{:3}: {}",
                    start + 1,
                    seq.substring(pos, end)
                )?;
            } else {
                writeln!(output, "{}", seq.substring(pos, end))?;
            }
        }
    }

    // Read the rest of the file into memory
    //let mut contents: Vec<u8> = vec![];
    //sa_file.read_to_end(&mut contents)?;
    //dbg!(&contents);

    // Determine how many bits to take
    //let mult = if (n as u64) < u32::MAX as u64 { 4 } else { 8 };
    //dbg!(&mult);

    //let raw_sa: Vec<_> = contents.drain(0..(n * mult)).collect();
    //dbg!(&raw_sa[..20]);

    //let chunked: Vec<_> = raw_sa.chunks(mult).collect();
    //dbg!(&chunked[..20]);

    //let converted: Vec<_> = chunked
    //    .iter()
    //    .map(|&v| usize::from_ne_bytes(v.into()))
    //    .collect();

    //.chunk_by(|chunk| usize::from_ne_bytes(chunk))
    //.collect();

    ////dbg!(&seq);

    //for pos in suffix_array {
    //    let p = pos as usize;
    //    println!("{pos} = {}", seq[p]);
    //    break;
    //}

    //let mut lcp: Vec<u8> = vec![0; *n as usize];
    //file.read_exact(&mut lcp)?;
    //dbg!(&lcp);

    Ok(())
}

fn read_suffix_array(filename: &str) -> Result<Vec<usize>> {
    let mut sa_file =
        File::open(&filename).map_err(|e| anyhow!("{filename}: {e}"))?;

    // The first 64-bits of the file contain the size of the SA
    let mut buffer = [0; 8];
    sa_file.read_exact(&mut buffer)?;

    // Convert the Vec<u8> to a usize
    let sa_len = usize::from_ne_bytes(buffer);

    // Allocate a buffer to hold the data
    //let size = if sa_len < u32::MAX as usize {
    //    mem::size_of::<u32>()
    //} else {
    //    mem::size_of::<u64>()
    //};
    //let mut buffer = vec![0u8; sa_len * size];

    // Read the bytes into the buffer
    let mut buffer = vec![0u8; sa_len * mem::size_of::<usize>()];
    sa_file.read_exact(&mut buffer)?;

    // How can I have either 32 or 64 vectors?
    // Convert the buffer into a slice of i32 integers
    //let suffix_array: &[u32] = unsafe {
    //    std::slice::from_raw_parts(buffer.as_ptr() as *const u32, sa_len)
    //};

    let suffix_array: &[usize] = unsafe {
        std::slice::from_raw_parts(buffer.as_ptr() as *const usize, sa_len)
    };

    Ok(suffix_array.to_vec())
}
