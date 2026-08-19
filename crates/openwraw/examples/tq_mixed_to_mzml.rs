//! Convert a mixed Waters TQ acquisition to mzML.
//!
//! Broad Q3 scans are written as spectra; MRM channels are written as
//! canonical SRM chromatograms with precursor/product m/z metadata.
//!
//! Usage:
//!
//! ```text
//! cargo run -p openwraw --example tq_mixed_to_mzml --release -- \
//!   path/to/bundle.raw out.mzML native-ms2
//!
//! cargo run -p openwraw --example tq_mixed_to_mzml --release -- \
//!   path/to/bundle.raw out.mzML pseudo-ms1
//! ```
//!
//! Add `--indexed` to emit indexed mzML.

use std::env;
use std::fs::File;
use std::io::BufWriter;
use std::time::Instant;

use openwraw::raw::tq_mrm::TqMrmReader;
use openwraw::raw::tq_reader::TqReader;
use openwraw::tq_mixed_mzml::{write_tq_mixed_indexed_mzml, write_tq_mixed_mzml};
use openwraw::tq_mzml::TqQ3MzmlMode;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = env::args().collect();
    if args.len() < 4 {
        eprintln!(
            "usage: tq_mixed_to_mzml <bundle.raw> <out.mzML> <native-ms2|pseudo-ms1> [--indexed]"
        );
        std::process::exit(2);
    }

    let bundle = &args[1];
    let output = &args[2];
    let mode = match args[3].as_str() {
        "native-ms2" => TqQ3MzmlMode::NativeMs2,
        "pseudo-ms1" => TqQ3MzmlMode::PseudoMs1,
        other => {
            eprintln!("unsupported mode {other:?}; use native-ms2 or pseudo-ms1");
            std::process::exit(2);
        }
    };
    let indexed = args.iter().any(|arg| arg == "--indexed");

    let q3 = TqReader::open(bundle)?;
    let mrm = TqMrmReader::open(bundle)?;
    eprintln!(
        "Detected {} broad Q3 function(s), {} Q3 scan(s), {} MRM function(s), {} MRM transition(s)",
        q3.q3_functions.len(),
        q3.total_q3_scan_count(),
        mrm.functions.len(),
        mrm.transition_count()
    );

    if q3.q3_functions.is_empty() && mrm.functions.is_empty() {
        return Err("no supported TQ Q3 or MRM functions found in bundle".into());
    }

    let started = Instant::now();
    let file = File::create(output)?;
    let mut writer = BufWriter::new(file);
    if indexed {
        write_tq_mixed_indexed_mzml(bundle, &mut writer, mode)?;
    } else {
        write_tq_mixed_mzml(bundle, &mut writer, mode)?;
    }

    let q3_label = match mode {
        TqQ3MzmlMode::NativeMs2 => "Q3=native-MS2",
        TqQ3MzmlMode::PseudoMs1 => "Q3=pseudo-MS1",
    };
    let indexed_label = if indexed { " indexed" } else { "" };
    eprintln!(
        "Wrote{indexed_label} mzML to {output} in {:.1}s ({q3_label}, MRM=SRM chromatograms)",
        started.elapsed().as_secs_f64()
    );
    Ok(())
}
