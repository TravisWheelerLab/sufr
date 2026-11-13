use anyhow::Result;
use clap::Parser;
use log::info;
use std::{
    fs::File,
    io::{self, BufWriter},
};
use sufr::{Cli, Command, LogLevel};

fn main() {
    if let Err(e) = run(Cli::parse()) {
        // the reader of stdout closed it, as head does once it has its
        // lines; what was written is correct, so this is not a failure
        let broken_pipe = e
            .chain()
            .filter_map(|cause| cause.downcast_ref::<io::Error>())
            .any(|e| e.kind() == io::ErrorKind::BrokenPipe);
        if broken_pipe {
            return;
        }
        eprintln!("Error: {e}");
        std::process::exit(1);
    }
}

fn run(cli: Cli) -> Result<()> {
    let Some(command) = cli.command else {
        unreachable!()
    };
    let args = command.run_args();
    env_logger::Builder::new()
        .filter_level(match args.log {
            Some(LogLevel::Debug) => log::LevelFilter::Debug,
            Some(LogLevel::Info) => log::LevelFilter::Info,
            _ => log::LevelFilter::Off,
        })
        .target(match args.log_file {
            Some(ref filename) => env_logger::Target::Pipe(Box::new(BufWriter::new(
                File::create(filename)?,
            ))),
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

    match &command {
        Command::Count(args) => sufr::count(args),
        Command::Create(args) => sufr::create(args),
        Command::Extract(args) => sufr::extract(args),
        Command::List(args) => sufr::list(args),
        Command::Locate(args) => sufr::locate(args),
        Command::Summarize(args) => sufr::summarize(args),
    }
}
