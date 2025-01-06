use anyhow::Result;
use clap::Parser;
use log::info;
use sufr::{Cli, Command, LogLevel};
use std::{io::BufWriter, fs::File};

// --------------------------------------------------
fn main() {
    if let Err(e) = run(Cli::parse()) {
        eprintln!("Error: {e}");
        std::process::exit(1);
    }
}

// --------------------------------------------------
fn run(args: Cli) -> Result<()> {
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

    match &args.command {
        // TODO: enable once "check" is usable
        //Some(Command::Check(args)) => {
        //    sufr::check(args)?;
        //    Ok(())
        //}
        Some(Command::Count(args)) => {
            sufr::count(args)?;
            Ok(())
        }
        Some(Command::Create(args)) => {
            sufr::create(args)?;
            Ok(())
        }
        Some(Command::Extract(args)) => {
            sufr::extract(args)?;
            Ok(())
        }
        Some(Command::List(args)) => {
            sufr::list(args)?;
            Ok(())
        }
        Some(Command::Locate(args)) => {
            sufr::locate(args)?;
            Ok(())
        }
        Some(Command::Summarize(args)) => {
            sufr::summarize(args)?;
            Ok(())
        }
        _ => unreachable!(),
    }
}
