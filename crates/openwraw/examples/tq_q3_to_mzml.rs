//! Convert broad Q3 scans from a Waters triple-quadrupole `.raw/` bundle to mzML.
//!
//! This example intentionally does not export MRM functions yet. It supports
//! two explicit semantic projections:
//!
//! ```text
//! cargo run -p openwraw --example tq_q3_to_mzml --release -- \
//!   path/to/bundle.raw out.mzML native-ms2
//!
//! cargo run -p openwraw --example tq_q3_to_mzml --release -- \
//!   path/to/bundle.raw out.mzML pseudo-ms1
//! ```
//!
//! Add `--indexed` to emit indexed mzML.

use std::env;
use std::fs::File;
use std::io::BufWriter;
use std::time::Instant;

use openwraw::tq_mzml::{write_tq_q3_indexed_mzml, write_tq_q3_mzml, TqQ3MzmlMode};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = env::args().collect();
    if args.len() < 4 {
        eprintln!(
            "usage: tq_q3_to_mzml <bundle.raw> <out.mzML> <native-ms2|pseudo-ms1> [--indexed]"
        );
        std::process::exit(2);
    }

    let bundle = &args[1];
    let out_path = &args[2];
    let mode = match args[3].as_str() {
        "native-ms2" => TqQ3MzmlMode::NativeMs2,
        "pseudo-ms1" => TqQ3MzmlMode::PseudoMs1,
        other => {
            eprintln!("unsupported mode {other:?}; use native-ms2 or pseudo-ms1");
            std::process::exit(2);
        }
    };
    let indexed = args.iter().any(|arg| arg == "--indexed");

    let reader = openwraw::raw::tq_reader::TqReader::open(bundle)?;
    eprintln!(
        "Detected {} broad Q3 function(s), {} Q3 scan(s), {} MRM function(s) not exported by this path",
        reader.q3_functions.len(),
        reader.total_q3_scan_count(),
        reader.mrm_function_count
    );

    if reader.q3_functions.is_empty() {
        return Err("no supported broad Q3 functions found in bundle".into());
    }

    let started = Instant::now();
    let file = File::create(out_path)?;
    let mut writer = BufWriter::new(file);
    if indexed {
        write_tq_q3_indexed_mzml(bundle, &mut writer, mode)?;
    } else {
        write_tq_q3_mzml(bundle, &mut writer, mode)?;
    }

    let semantic_label = match mode {
        TqQ3MzmlMode::NativeMs2 => "native MS2/Q3 semantics",
        TqQ3MzmlMode::PseudoMs1 => "explicit pseudo-MS1 projection",
    };
    let indexed_label = if indexed { " indexed" } else { "" };
    eprintln!(
        "Wrote{indexed_label} mzML to {out_path} in {:.1}s ({semantic_label})",
        started.elapsed().as_secs_f64()
    );
    Ok(())
}
