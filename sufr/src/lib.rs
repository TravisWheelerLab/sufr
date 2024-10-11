use anyhow::{anyhow, bail, Result};
use clap::{builder::PossibleValue, Parser, ValueEnum};
use format_num::NumberFormat;
use libsufr::{
    read_sequence_file, read_suffix_length, FromUsize, Int, SufrBuilder, SufrFile,
};
use log::info;
use regex::Regex;
use std::{
    cmp::min,
    ffi::OsStr,
    fmt::Debug,
    fs::File,
    io::{self, BufWriter, Write},
    num::NonZeroUsize,
    ops::Range,
    path::PathBuf,
    time::Instant,
};
//use u4::{AsNibbles, U4x2, U4};

// --------------------------------------------------
#[derive(Parser, Debug)]
#[command(arg_required_else_help = true, version, about)]
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
        let sufr_file: SufrFile<u32> = SufrFile::read(&args.filename)?;
        _check(sufr_file, args)
    } else {
        let sufr_file: SufrFile<u64> = SufrFile::read(&args.filename)?;
        _check(sufr_file, args)
    }
}

// --------------------------------------------------
pub fn _check<T>(sufr_file: SufrFile<T>, _args: &CheckArgs) -> Result<()>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + Debug,
{
    //dbg!(&sufr_file);
    //println!("text         = {}", sufr_file.get_text()?);
    //println!("raw text     = {:?}", sufr_file.get_raw_text()?);
    //println!("suffix_array = {:?}", sufr_file.get_suffix_array()?);
    //println!("lcp          = {:?}", sufr_file.get_lcp()?);

    //let suffix_array_iter = sufr_file.get_suffix_array_iter();
    //for i in 0..sufr_file.num_suffixes.to_usize() {
    //    let val = suffix_array_iter.get(i)?;
    //    println!("SA {i:3}: {val}");
    //}
    ////let r = &suffix_array_iter.get_range(1..5)?;
    ////dbg!(r);
    //
    //for (i, val) in suffix_array_iter.enumerate() {
    //    println!("suffix_array {i:3}: {val}");
    //}
    let now = Instant::now();
    let errors = sufr_file.check_suffix_array()?;
    let num_errors = errors.len();
    let num_fmt = NumberFormat::new();
    println!(
        "Checked {} suffix{}, found {} error{} in suffix array.",
        num_fmt.format(",.0", sufr_file.num_suffixes.to_usize() as f64),
        if sufr_file.num_suffixes.to_usize() == 1 {
            ""
        } else {
            "es"
        },
        num_fmt.format(",.0", num_errors as f64),
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
        let sa: SufrBuilder<u32> = SufrBuilder::new(
            seq_data.seq,
            len as u32,
            args.max_context.map(|val| val as u32),
            args.is_dna,
            start_positions,
            seq_data.headers,
            args.num_partitions,
        )?;
        _create(sa, args, now)
    } else {
        let start_positions: Vec<u64> =
            seq_data.start_positions.iter().map(|&v| v as u64).collect();
        let sa: SufrBuilder<u64> = SufrBuilder::new(
            seq_data.seq,
            len,
            args.max_context.map(|val| val as u64),
            args.is_dna,
            start_positions,
            seq_data.headers,
            args.num_partitions,
        )?;
        _create(sa, args, now)
    }
}

// --------------------------------------------------
// Helper for "create" that checks/writes SA/LCP
pub fn _create<T>(sa: SufrBuilder<T>, args: &CreateArgs, timer: Instant) -> Result<()>
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
    let mut output =
        BufWriter::new(File::create(outfile).map_err(|e| anyhow!("{outfile}: {e}"))?);
    let now = Instant::now();
    let bytes_written = sa.write(&mut output)?;
    let num_fmt = NumberFormat::new();
    info!(
        "Wrote {} byte{} to '{outfile}' in {:?}",
        num_fmt.format(",.0", bytes_written as f64),
        if bytes_written == 1 { "" } else { "s" },
        now.elapsed()
    );
    info!("Total time: {:?}", timer.elapsed());

    Ok(())
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
pub fn extract(args: &ExtractArgs) -> Result<()> {
    let len = read_suffix_length(&args.filename)? as u64;
    if len < u32::MAX as u64 {
        let now = Instant::now();
        let sufr_file: SufrFile<u32> = SufrFile::read(&args.filename)?;
        info!("Read sufr file in {:?}", now.elapsed());
        _extract(sufr_file, args)
    } else {
        let now = Instant::now();
        let sufr_file: SufrFile<u64> = SufrFile::read(&args.filename)?;
        info!("Read sufr file in {:?}", now.elapsed());
        _extract(sufr_file, args)
    }
}

// --------------------------------------------------
pub fn _extract<T>(sufr_file: SufrFile<T>, args: &ExtractArgs) -> Result<()>
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

    let text_len = sufr_file.len.to_usize();
    let width = text_len.to_string().len() + 1;
    for pos in positions {
        if let Some(suffix) = sufr_file.suffix_array.get(pos).map(|v| v.to_usize()) {
            let end = min(args.max_len.map_or(text_len, |len| suffix + len), text_len);
            let seq = String::from_utf8(sufr_file.text[suffix..end].to_vec())?;

            if args.number {
                writeln!(output, "{suffix:width$}: {seq}")?;
            } else {
                writeln!(output, "{seq}")?;
            }
        } else {
            break;
        }
    }
    info!("Extracted suffixes in {:?}", now.elapsed());

    Ok(())
}

// --------------------------------------------------
pub fn search(args: &SearchArgs) -> Result<()> {
    let len = read_suffix_length(&args.filename)? as u64;
    if len < u32::MAX as u64 {
        let now = Instant::now();
        let sa: SufrFile<u32> = SufrFile::read(&args.filename)?;
        info!("Read sufr file in {:?}", now.elapsed());
        _search(sa, args)
    } else {
        let now = Instant::now();
        let sa: SufrFile<u64> = SufrFile::read(&args.filename)?;
        info!("Read sufr file in {:?}", now.elapsed());
        _search(sa, args)
    }
}

// --------------------------------------------------
pub fn _search<T>(sufr_file: SufrFile<T>, args: &SearchArgs) -> Result<()>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + Debug,
{
    let now = Instant::now();
    let res = match sufr_file.search(&args.query) {
        Some((start, stop)) => format!("found in range {}..{}", start + 1, stop + 1),
        _ => "not found".to_string(),
    };
    println!("Query '{}' {res} in {:?}", args.query, now.elapsed());
    Ok(())
}
