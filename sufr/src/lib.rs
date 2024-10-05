use anyhow::{anyhow, bail, Result};
use clap::{builder::PossibleValue, Parser, ValueEnum};
use format_num::NumberFormat;
use libsufr::{
    read_sequence_file, read_suffix_length, FromUsize, Int, SuffixArray,
    SuffixArrayBuilder,
};
use log::{debug, info};
use regex::Regex;
use std::{
    ffi::OsStr,
    fmt::Debug,
    fs::File,
    io::{self, BufWriter, Write},
    ops::Range,
    path::PathBuf,
    time::Instant,
};
//use u4::{AsNibbles, U4x2, U4};

// --------------------------------------------------
#[derive(Parser, Debug)]
#[command(arg_required_else_help = true)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Option<Command>,

    /// Log level
    #[arg(short, long)]
    pub log: Option<LogLevel>,

    /// Log file
    #[arg(long)]
    pub log_file: Option<String>,
}

#[derive(Parser, Debug)]
pub enum Command {
    /// Create suffix array
    Create(CreateArgs),

    /// Check correctness of suffix array/LCP
    Check(CheckArgs),

    /// Read suffix array and extract sequences
    Extract(ExtractArgs),

    /// Search a suffix array
    Search(SearchArgs),
}

#[derive(Debug, Parser)]
#[command(author, version, about)]
pub struct CheckArgs {
    /// Sufr file
    #[arg(value_name = "SUFR")]
    pub filename: String,

    /// List errors
    #[arg(short, long)]
    pub verbose: bool,
}

#[derive(Debug, Parser)]
#[command(author, version, about)]
pub struct CreateArgs {
    /// Input file
    #[arg(value_name = "INPUT")]
    pub input: String,

    /// Subproblem count
    #[arg(short, long, value_name = "NUM_PARTS", default_value = "16")]
    pub num_partitions: usize,

    /// Max context
    #[arg(short, long, value_name = "CONTEXT")]
    pub max_context: Option<usize>,

    /// Number of threads
    #[arg(short, long, value_name = "THREADS")]
    pub threads: Option<usize>,

    /// Output file
    #[arg(short, long, value_name = "OUTPUT")]
    pub output: Option<String>,

    /// Input is DNA, ignore sequences starting with 'N'
    #[arg(short('d'), long("dna"))]
    pub is_dna: bool,
}

/// Extract suffixes from sufr file
#[derive(Debug, Parser)]
#[command(author, version, about)]
pub struct ExtractArgs {
    /// Sufr file
    #[arg(value_name = "SUFR")]
    pub filename: String,

    /// Maximum length of sequence
    #[arg(short, long, value_name = "MAX")]
    pub max_len: Option<usize>,

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

/// Search a suffix array
#[derive(Debug, Parser)]
#[command(author, version, about)]
pub struct SearchArgs {
    /// Query
    #[arg(value_name = "QUERY")]
    pub query: String,

    /// Sufr file
    #[arg(value_name = "SUFR")]
    pub filename: String,
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
    let len = read_suffix_length(&args.filename)? as u64;
    if len < u32::MAX as u64 {
        let sa: SuffixArray<u32> = SuffixArray::read(&args.filename)?;
        _check(sa, args)
    } else {
        let sa: SuffixArray<u64> = SuffixArray::read(&args.filename)?;
        _check(sa, args)
    }
}

// --------------------------------------------------
pub fn _check<T>(sa: SuffixArray<T>, _args: &CheckArgs) -> Result<()>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + Debug,
{
    let now = Instant::now();

    // Manual for now
    let num_errors = 0;
    if sa.suffix_array.len() > 1 {
        for i in 1..sa.suffix_array.len() {
            let cur = sa.string_at(sa.suffix_array[i].to_usize());
            let prev = sa.string_at(sa.suffix_array[i-1].to_usize());
            if prev > cur {
                eprintln!("{i}: '{cur}' vs '{prev}'");
            }
        }
    }

    println!(
        "Found {num_errors} error{} in suffix array.",
        if num_errors == 1 { "" } else { "s" },
    );

    //let errors = sa.check_order();
    //let num_errors = errors.len();
    //
    //if args.verbose {
    //    for (i, pos) in errors.iter().enumerate() {
    //        println!(
    //            "{:3}: pos {pos:5} {}",
    //            i + 1,
    //            sa.string_at(pos.to_usize())
    //        );
    //    }
    //}
    //
    //let lcp_errors = sa.check_lcp();
    //let num_errors = lcp_errors.len();
    //println!(
    //    "Found {num_errors} error{} in LCP",
    //    if num_errors == 1 { "" } else { "s" },
    //);
    //
    //if args.verbose {
    //    for (i, pos) in lcp_errors.iter().enumerate() {
    //        println!("{:3}: pos {pos:5}", i + 1,);
    //    }
    //}
    //
    println!("Finished checking in {:?}.", now.elapsed());
    Ok(())
}

// --------------------------------------------------
pub fn create(args: &CreateArgs) -> Result<()> {
    let num_threads = args.threads.unwrap_or(num_cpus::get());
    info!(
        "Using {num_threads} thread{}",
        if num_threads == 1 { "" } else { "s" }
    );
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .unwrap();

    // Read sequence input
    let now = Instant::now();
    let seq_data = read_sequence_file(&args.input)?;
    let len = seq_data.seq.len() as u64;
    let num_fmt = NumberFormat::new();
    info!(
        "Read raw input of len {} in {:?}",
        num_fmt.format(",.0", len as f64),
        now.elapsed()
    );
    //debug!("Raw input '{:?}'", seq_data.seq);

    if len < u32::MAX as u64 {
        let start_positions: Vec<u32> =
            seq_data.start_positions.iter().map(|&v| v as u32).collect();
        let sa: SuffixArrayBuilder<u32> = SuffixArrayBuilder::new(
            seq_data.seq,
            len as u32,
            args.max_context.map(|val| val as u32),
            args.is_dna,
            start_positions,
            seq_data.headers,
            args.num_partitions,
        )?;
        _create(sa, args)
    } else {
        let start_positions: Vec<u64> =
            seq_data.start_positions.iter().map(|&v| v as u64).collect();
        let sa: SuffixArrayBuilder<u64> = SuffixArrayBuilder::new(
            seq_data.seq,
            len,
            args.max_context.map(|val| val as u64),
            args.is_dna,
            start_positions,
            seq_data.headers,
            args.num_partitions,
        )?;
        _create(sa, args)
    }
}

// --------------------------------------------------
// Helper for "create" that checks/writes SA/LCP
pub fn _create<T>(sa: SuffixArrayBuilder<T>, args: &CreateArgs) -> Result<()>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + Debug,
{
    let outfile = &args.output.clone().unwrap_or(format!(
        "{}.sufr",
        PathBuf::from(&args.input)
            .file_stem()
            .unwrap_or(OsStr::new("out"))
            .to_string_lossy()
    ));
    let mut output = BufWriter::new(
        File::create(outfile).map_err(|e| anyhow!("{outfile}: {e}"))?,
    );
    let now = Instant::now();
    let bytes_written = sa.write(&mut output)?;
    let num_fmt = NumberFormat::new();
    info!(
        "Wrote {} byte{} to '{outfile}' in {:?}",
        num_fmt.format(",.0", bytes_written as f64),
        if bytes_written == 1 { "" } else { "s" },
        now.elapsed()
    );

    Ok(())
}

// --------------------------------------------------
// Parse an index from a string representation of an integer.
// Ensures the number does not start with '+'.
fn parse_index(input: &str) -> Result<usize> {
    let value_error = || anyhow!(r#"illegal list value: "{input}""#);
    input
        .starts_with('+')
        .then(|| Err(value_error()))
        .unwrap_or_else(|| input.parse::<usize>().map_err(|_| value_error()))
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
pub fn extract(args: &ExtractArgs) -> Result<()> {
    let len = read_suffix_length(&args.filename)? as u64;
    if len < u32::MAX as u64 {
        let now = Instant::now();
        let sa: SuffixArray<u32> = SuffixArray::read(&args.filename)?;
        info!("Read sufr file in {:?}", now.elapsed());
        _extract(sa, args)
    } else {
        let now = Instant::now();
        let sa: SuffixArray<u64> = SuffixArray::read(&args.filename)?;
        info!("Read sufr file in {:?}", now.elapsed());
        _extract(sa, args)
    }
}

// --------------------------------------------------
pub fn _extract<T>(sa: SuffixArray<T>, args: &ExtractArgs) -> Result<()>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + Debug,
{
    let now = Instant::now();
    let positions: Vec<_> = parse_pos(&args.extract)?
        .into_iter()
        .flat_map(|r| r.collect::<Vec<_>>())
        .collect();

    let mut output: Box<dyn Write> = match &args.output {
        Some(out_name) => Box::new(File::create(out_name)?),
        _ => Box::new(io::stdout()),
    };

    let sa_len = sa.text.len();
    for start in positions {
        if let Some(pos) = sa.suffix_array.get(start).map(|v| v.to_usize()) {
            let mut end = args.max_len.map_or(sa_len, |len| pos + len);
            if end > sa_len {
                end = sa_len
            }
            let seq = String::from_utf8(sa.text[pos..end].to_vec())?;
            if args.number {
                writeln!(output, "{pos:3}: {seq}")?;
            } else {
                writeln!(output, "{seq}")?;
            }
        }
    }
    info!("Extracted suffixes in {:?}", now.elapsed());

    Ok(())
}

//// --------------------------------------------------
//pub fn search(args: &SearchArgs) -> Result<()> {
//    let len = read_suffix_length(&args.filename)? as u64;
//    if len < u32::MAX as u64 {
//        let now = Instant::now();
//        let sa: SuffixArray<u32> = SuffixArray::read(&args.filename)?;
//        info!("Read sufr file in {:?}", now.elapsed());
//        _search(sa, args)
//    } else {
//        let now = Instant::now();
//        let sa: SuffixArray<u64> = SuffixArray::read(&args.filename)?;
//        info!("Read sufr file in {:?}", now.elapsed());
//        _search(sa, args)
//    }
//}
//
//// --------------------------------------------------
//pub fn _search<T>(sa: SuffixArrayBuilder<T>, args: &SearchArgs) -> Result<()>
//where
//    T: Int + FromUsize<T> + Sized + Send + Sync + Debug,
//{
//    let now = Instant::now();
//    let res = match sa.search(&args.query) {
//        Some(range) => format!("found in range {range:?}"),
//        _ => "not found".to_string(),
//    };
//    println!("Query '{}' {res} in {:?}", args.query, now.elapsed());
//    Ok(())
//}
