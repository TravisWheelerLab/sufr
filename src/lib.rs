mod suffix_array;

use anyhow::{anyhow, bail, Result};
use clap::{builder::PossibleValue, Parser, ValueEnum};
use log::info;
use regex::Regex;
use std::{
    fs::{self, File},
    io::{self, BufRead, BufReader, BufWriter, Read, Write},
    mem,
    num::NonZeroUsize,
    ops::Range,
    path::PathBuf,
    //slice,
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
    #[clap(alias = "cr")]
    Create(CreateArgs),

    /// Read suffix array and extract sequences
    #[clap(alias = "re")]
    Read(ReadArgs),
}

#[derive(Debug, Parser)]
#[command(author, version, about)]
pub struct CreateArgs {
    /// Input file
    #[arg(value_name = "INPUT")]
    pub input: String,

    /// Search pattern
    #[arg(value_name = "PATTERN")]
    pub pattern: Option<String>,

    /// Subproblem count
    #[arg(short, long, value_name = "SUBPROBLEMS", default_value = "8192")]
    pub subproblem_count: usize,

    /// Max context
    #[arg(short, long, value_name = "CONTEXT")]
    pub max_context: Option<usize>,

    /// Number of threads
    #[arg(short, long, value_name = "THREADS")]
    pub threads: Option<usize>,

    /// Output file
    #[arg(short, long, value_name = "OUTPUT", default_value = "sufr.sa")]
    pub output: String,

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
    pub array: PathBuf,

    /// Sequence file
    #[arg(short, long, value_name = "SEQ")]
    pub sequence: PathBuf,

    /// Maximum length of sequence
    #[arg(short, long, value_name = "MAX", default_value = "0")]
    pub max_len: usize,

    /// Extract positions
    #[arg(short, long, value_name = "EXTRACT", default_value = "1")]
    pub extract: String,
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

    let start = Instant::now();
    let suf_arr = SuffixArray::new(open(&args.input)?, args.max_context);
    info!(
        "Finished reading input of len {} in {:?}",
        suf_arr.len,
        start.elapsed()
    );
    //debug!("Raw input '{:?}'", suf_arr.text);

    let start = Instant::now();
    //let sa = suf_arr.generate();
    let sufs = suf_arr.sort_subarrays(args.subproblem_count);
    dbg!(&sufs);
    let suffixes: Vec<_> = sufs.iter().map(|t| t.0.clone()).collect();
    let pivots = suf_arr.select_pivots(&suffixes);
    dbg!(&pivots);

    let pivot_locs = suf_arr.locate_pivots(&suffixes, &pivots);
    dbg!(&pivot_locs);

    let part_subs = suf_arr.partition_subarrays(&suffixes, &pivot_locs);
    dbg!(&part_subs);

    let merged = suf_arr.merge_part_subs(&part_subs);
    dbg!(merged);

    info!("Generated suffix array in {:?}", start.elapsed());
    //dbg!(&sa);

    let mut out = File::create(&args.output)?;

    // Dynamically choose 32/64
    //let size = if suf_arr.len < u32::MAX as usize {
    //    mem::size_of::<u32>()
    //} else {
    //    mem::size_of::<u64>()
    //};
    //let slice_u8: &[u8] = unsafe {
    //    slice::from_raw_parts(sa.as_ptr() as *const _, suf_arr.len * size)
    //};

    // Write out suffix array length
    out.write(&usize_to_bytes(suf_arr.len))?;

    // Write out suffix array as raw bytes
    //let slice_u8: &[u8] = unsafe {
    //    slice::from_raw_parts(
    //        sa.as_ptr() as *const _,
    //        suf_arr.len * mem::size_of::<usize>(),
    //    )
    //};
    //out.write_all(slice_u8)?;

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
    let mut sa_file = File::open(&args.array)?;

    // The first 64-bits of the file contain the size of the SA
    let mut buffer = [0; 8];
    sa_file.read_exact(&mut buffer)?;

    // Convert the Vec<u8> to a usize
    let sa_len = usize::from_ne_bytes(buffer);
    println!("Length of suffix array: {sa_len}");

    let start = Instant::now();

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

    println!("Read SA in {:?}", start.elapsed());

    let seq = fs::read_to_string(&args.sequence)?.to_uppercase();
    let seq = format!("{}$", seq.trim());

    let positions: Vec<_> = parse_pos(&args.extract)?
        .into_iter()
        .flat_map(|r| r.collect::<Vec<_>>())
        .collect();

    for start in positions {
        if let Some(&pos) = suffix_array.get(start) {
            //let pos = pos as usize;
            let end = if args.max_len > 0 {
                pos + args.max_len
            } else {
                sa_len
            };

            println!(
                "{:3}: {pos:3} = {}",
                start + 1,
                seq.substring(pos, end)
            );
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
