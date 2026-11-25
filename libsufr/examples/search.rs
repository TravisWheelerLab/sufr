use libsufr::SuffixArray;

type Result<T, E = Box<dyn std::error::Error>> = std::result::Result<T, E>;

fn main() -> Result<()> {
    // the text is ACGTNNACGT$
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../data/inputs/1.sufr");
    let suffix_array = SuffixArray::open(path)?;

    let hits = suffix_array.search(b"ACG")?;
    println!("count: {}", hits.count());
    println!("locate: {:?}", hits.locate()?);

    for extraction in hits.extract(1, Some(3))? {
        let start = extraction.location.sequence_start;
        let window = extraction.sequence_range;
        let text = suffix_array.text(start + window.start..start + window.end)?;
        println!(
            "extract: {} with the hit at {}",
            String::from_utf8_lossy(&text),
            extraction.suffix_offset
        );
    }

    for hits in suffix_array.search_all(["AC", "GT", "XX"])? {
        println!("count: {}", hits.count());
    }

    // the same queries, spread over the rayon pool
    let queries = vec![String::from("AC"), String::from("GT")];
    for hits in suffix_array.par_search_all(&queries)? {
        println!("count: {}", hits.count());
    }

    Ok(())
}
