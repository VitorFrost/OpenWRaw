# Triple-quadrupole / low-resolution Waters RAW

> Experimental. This path is intentionally separate from the established QTOF/IMS `Reader` while TQ support is validated across independently redistributable datasets.

OpenWRaw can represent mixed Waters triple-quadrupole acquisitions containing broad Q3 scans and MRM functions without requiring Waters SDK libraries.

## Data model

The TQ path keeps the two acquisition families distinct:

- **Broad Q3 scans** are decoded as dense m/z-intensity spectra from the low-resolution 6-byte packed representation.
- **MRM functions** are decoded from 4-byte packed intensities and emitted as one selected-reaction-monitoring chromatogram per Q1 -> Q3 transition.

The mixed mzML writer therefore produces:

```text
mzML
├── spectrumList
│   └── broad Q3 scans
└── chromatogramList
    └── MRM / SRM transition chromatograms
```

MRM chromatograms use PSI-MS `MS:1001473` (`selected reaction monitoring chromatogram`) and carry explicit precursor and product m/z values.

## Q3 semantics

Waters triple-quadrupole broad scans may be recorded as MS2/Q3 scans even when the resulting broad mass spectrum is intended for untargeted downstream processing. OpenWRaw therefore requires an explicit choice rather than silently changing the MS level.

### `native-ms2`

Preserves the vendor acquisition semantics. Broad Q3 scans remain MS2 and retain the function-level precursor/set-mass context.

Use this for archival conversion and format-faithful interchange.

### `pseudo-ms1`

Projects the broad Q3 spectrum to MS1 for downstream tools that require an MS1 survey stream. The transformation is explicit and the emitted spectrum filter metadata identifies it as an OpenWRaw pseudo-MS1 projection.

MRM data is unaffected by this choice and remains SRM chromatograms.

## Headless mixed conversion

From a source checkout:

```bash
cargo run -p openwraw --example tq_mixed_to_mzml --release -- \
  path/to/bundle.raw output.mzML native-ms2
```

For an untargeted workflow that needs the broad Q3 scans presented as MS1:

```bash
cargo run -p openwraw --example tq_mixed_to_mzml --release -- \
  path/to/bundle.raw output.mzML pseudo-ms1
```

Add `--indexed` to either command for indexed mzML:

```bash
cargo run -p openwraw --example tq_mixed_to_mzml --release -- \
  path/to/bundle.raw output.mzML pseudo-ms1 --indexed
```

The converter reports the number of broad-Q3 functions/scans and MRM functions/transitions it detected before writing output.

## Q3-only diagnostic conversion

A separate example excludes MRM chromatograms and is useful while debugging the broad-scan decoder:

```bash
cargo run -p openwraw --example tq_q3_to_mzml --release -- \
  path/to/bundle.raw q3-only.mzML native-ms2
```

or:

```bash
cargo run -p openwraw --example tq_q3_to_mzml --release -- \
  path/to/bundle.raw q3-only.mzML pseudo-ms1
```

## Rust API

Broad Q3 scans can be accessed independently:

```rust
use openwraw::raw::tq_reader::TqReader;

let reader = TqReader::open("sample.raw")?;
for scan in reader.iter_scans() {
    let scan = scan?;
    println!("function={} rt={} points={}",
        scan.function_index,
        scan.retention_time_min,
        scan.spectrum.mz.len());
}
```

MRM traces have a separate reader:

```rust
use openwraw::raw::tq_mrm::TqMrmReader;

let reader = TqMrmReader::open("sample.raw")?;
for trace in reader.chromatograms()? {
    println!("{} -> {}: {} points",
        trace.precursor_mz,
        trace.product_mz,
        trace.intensity.len());
}
```

For mixed mzML output:

```rust
use openwraw::tq_mixed_mzml::write_tq_mixed_mzml;
use openwraw::tq_mzml::TqQ3MzmlMode;

write_tq_mixed_mzml(
    "sample.raw",
    &mut output,
    TqQ3MzmlMode::PseudoMs1,
)?;
```

## Confidentiality and regression fixtures

The TQ implementation was developed using format-level observations that can be expressed independently of any one acquisition. Closed or confidential RAW data must **not** be committed as repository fixtures, copied into documentation, or reproduced as run-specific hexadecimal examples.

Repository tests for this path must use either:

1. synthetic binary records constructed by the tests themselves, or
2. independently redistributable public datasets whose licensing and provenance allow their use.

Private data may be used locally to verify that a generic decoder behaves correctly, but tests and documentation committed to the repository must not contain sample identifiers, private method details, real transition lists, chromatographic results, internal paths, serial numbers, or byte sequences copied from a closed run.

## Current limitations

The mixed converter currently focuses on the information required for structurally correct Q3 spectra and MRM chromatograms. Additional Waters side-file metadata such as compound labels and per-transition method parameters can be exposed later once their API and mzML representation are defined cleanly.

The legacy QTOF/IMS `Reader::open` path is intentionally unchanged. TQ support currently uses the dedicated TQ readers and conversion functions documented above.
