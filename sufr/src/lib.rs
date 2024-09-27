use anyhow::{anyhow, bail, Result};
use clap::{builder::PossibleValue, Parser, ValueEnum};
use format_num::NumberFormat;
use libsufr::{
    read_sequence_file, read_suffix_length, FromUsize, Int, SuffixArray,
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

    /// Log level
    #[arg(short, long)]
    pub log: Option<LogLevel>,

    /// Log file
    #[arg(long)]
    pub log_file: Option<String>,
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

//#[derive(Debug)]
//pub struct PackedSeq {
//    data: Vec<u8>,
//}

//impl PackedSeq {
//    pub fn new(text: Vec<u8>) -> Self {
//        let u4s: Vec<_> = text.iter().filter_map(Self::byte_to_u4).collect();
//        let mut data: Vec<u8> = vec![0u8; text.len() / 2];
//        AsNibbles(&mut data).pack_from_slice(&u4s);
//        Self { data }
//    }

//    pub fn values(&self) -> AsNibbles<&Vec<u8>> {
//        AsNibbles(&self.data)
//    }

//    pub fn get(&self, index: usize) -> U4 {
//        let idx = index / 2;
//        let rem = index % 2;
//        let val = U4x2::from_byte(self.data[idx]);
//        if rem == 0 {
//            val.left()
//        } else {
//            val.right()
//        }
//    }

//    fn byte_to_u4(byte: &u8) -> Option<U4> {
//        match byte {
//            36 => Some(U4::Dec00), // $
//            65 => Some(U4::Dec01), // A
//            66 => Some(U4::Dec02), // B
//            67 => Some(U4::Dec03), // C
//            68 => Some(U4::Dec04), // D
//            71 => Some(U4::Dec05), // G
//            72 => Some(U4::Dec06), // H
//            75 => Some(U4::Dec07), // K
//            77 => Some(U4::Dec08), // M
//            78 => Some(U4::Dec09), // N
//            82 => Some(U4::Dec10), // R
//            83 => Some(U4::Dec11), // S
//            84 => Some(U4::Dec12), // T
//            86 => Some(U4::Dec13), // V
//            87 => Some(U4::Dec14), // W
//            89 => Some(U4::Dec15), // Y
//            _ => None,             // Anything else
//        }
//    }

//    fn u4_to_byte(val: U4) -> u8 {
//        match val {
//            U4::Dec00 => 36,
//            U4::Dec01 => 65,
//            U4::Dec02 => 66,
//            U4::Dec03 => 67,
//            U4::Dec04 => 68,
//            U4::Dec05 => 71,
//            U4::Dec06 => 72,
//            U4::Dec07 => 75,
//            U4::Dec08 => 77,
//            U4::Dec09 => 78,
//            U4::Dec10 => 82,
//            U4::Dec11 => 83,
//            U4::Dec12 => 84,
//            U4::Dec13 => 86,
//            U4::Dec14 => 87,
//            U4::Dec15 => 89,
//        }
//    }
//}

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
pub fn _check<T>(sa: SuffixArray<T>, args: &CheckArgs) -> Result<()>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + Debug,
{
    let now = Instant::now();
    let errors = sa.check_order();
    let num_errors = errors.len();

    if args.verbose {
        for (i, pos) in errors.iter().enumerate() {
            println!(
                "{:3}: pos {pos:5} {}",
                i + 1,
                sa.string_at(pos.to_usize())
            );
        }
    }

    println!(
        "Found {num_errors} error{} in suffix array.",
        if num_errors == 1 { "" } else { "s" },
    );

    let lcp_errors = sa.check_lcp();
    let num_errors = lcp_errors.len();
    println!(
        "Found {num_errors} error{} in LCP",
        if num_errors == 1 { "" } else { "s" },
    );

    if args.verbose {
        for (i, pos) in lcp_errors.iter().enumerate() {
            println!("{:3}: pos {pos:5}", i + 1,);
        }
    }

    println!("Finished checking in {:?}.", now.elapsed());
    Ok(())
}

// --------------------------------------------------
//pub fn read_input_u4(input: impl BufRead) -> PackedSeq {
//    //text.push(b'$');
//    let mut text: Vec<u8> = input
//        .bytes()
//        .map_while(Result::ok)
//        .filter(|&b| (0b1000001..=0b1011010).contains(&b)) // A-Z
//        .map(|b| b & 0b1011111) // Uppercase (mask w/32)
//        .collect();
//    text.push(b'$');

//    // The len of text must be even b/c we pack 2 u4's into each each u8
//    if text.len() % 2 == 1 {
//        text.push(b'$');
//    }
//    PackedSeq::new(text)
//}

// --------------------------------------------------
pub fn create(args: &CreateArgs) -> Result<()> {
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
    debug!("Raw input '{:?}'", seq_data.seq);

    //println!("{:?}", String::from_utf8(text.clone()));
    //println!("{}", headers.join(", "));
    //let packed = read_input_u4(&mut BufReader::new(File::open(&args.input)?));
    //println!("Input as U4 {:?}", &packed);
    //for (i, val) in packed.values().iter().enumerate() {
    //    println!(
    //        "{i}: {val:04b} => {} => {}",
    //        PackedSeq::u4_to_byte(val),
    //        packed.get(i),
    //    );
    //}

    if len < u32::MAX as u64 {
        let start_positions: Vec<u32> =
            seq_data.start_positions.iter().map(|&v| v as u32).collect();
        let sa: SuffixArray<u32> = SuffixArray::new(
            seq_data.seq,
            len as u32,
            args.max_context.map(|val| val as u32),
            args.is_dna,
            start_positions,
            seq_data.headers,
            args.num_partitions,
        );
        _create(sa, args)?;
    } else {
        let start_positions: Vec<u64> =
            seq_data.start_positions.iter().map(|&v| v as u64).collect();
        let sa: SuffixArray<u64> = SuffixArray::new(
            seq_data.seq,
            len,
            args.max_context.map(|val| val as u64),
            args.is_dna,
            start_positions,
            seq_data.headers,
            args.num_partitions,
        );
        _create(sa, args)?;
    };

    Ok(())
}

// --------------------------------------------------
// Helper for "create" that checks/writes SA/LCP
pub fn _create<T>(sa: SuffixArray<T>, args: &CreateArgs) -> Result<()>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + Debug,
{
    // Check
    if args.check {
        debug!("Sorted = {:?}", &sa.suffix_array);
        debug!(
            "Suffixes = {:#?}",
            sa.suffix_array
                .iter()
                .map(|v| sa.string_at(v.to_usize()))
                .collect::<Vec<_>>()
        );
        debug!("LCP = {:#?}", sa.lcp);
        let now = Instant::now();
        let num_errors = sa.check_order().len();
        info!(
            "Checked order, found {num_errors} error{} in {:?}",
            if num_errors == 1 { "" } else { "s" },
            now.elapsed()
        );
        let lcp_errors = sa.check_lcp();
        let num_errors = lcp_errors.len();
        info!(
            "Checked LCP, found {num_errors} error{} in {:?}",
            if num_errors == 1 { "" } else { "s" },
            now.elapsed()
        );
        if num_errors > 0 {
            debug!("{lcp_errors:?}");
        }
    }

    let now = Instant::now();
    let outfile = &args.output.clone().unwrap_or(format!(
        "{}.sufr",
        PathBuf::from(&args.input)
            .file_stem()
            .unwrap_or(OsStr::new("out"))
            .to_string_lossy()
    ));
    let out = BufWriter::new(
        File::create(outfile).map_err(|e| anyhow!("{outfile}: {e}"))?,
    );
    let bytes_written = sa.write(out)?;
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
        let sa: SuffixArray<u32> = SuffixArray::read(&args.filename)?;
        _extract(sa, args)
    } else {
        let sa: SuffixArray<u64> = SuffixArray::read(&args.filename)?;
        _extract(sa, args)
    }
}

// --------------------------------------------------
pub fn _extract<T>(sa: SuffixArray<T>, args: &ExtractArgs) -> Result<()>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + Debug,
{
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

    Ok(())
}

// --------------------------------------------------
pub fn search(args: &SearchArgs) -> Result<()> {
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

    let len = read_suffix_length(&args.filename)? as u64;
    if len < u32::MAX as u64 {
        let sa: SuffixArray<u32> = SuffixArray::read(&args.filename)?;
        _search(sa, args)
    } else {
        let sa: SuffixArray<u64> = SuffixArray::read(&args.filename)?;
        _search(sa, args)
    }
}

// --------------------------------------------------
pub fn _search<T>(sa: SuffixArray<T>, args: &SearchArgs) -> Result<()>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + Debug,
{
    let now = Instant::now();
    let res = sa.search(&args.query);
    println!("Found '{}' at {res:?} in {:?}", args.query, now.elapsed());
    Ok(())
}
