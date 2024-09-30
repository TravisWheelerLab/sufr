use anyhow::Result;
use clap::Parser;
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

    match &args.command {
        Some(Command::Check(args)) => {
            sufr::check(args)?;
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
        Some(Command::Search(args)) => {
            sufr::search(args)?;
            Ok(())
        }
        _ => unreachable!(),
    }
}
