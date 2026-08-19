# Triple-quadrupole / low-resolution Waters RAW

> Experimental. This path is intentionally separate from the established QTOF/IMS `Reader` while TQ support is validated across independently redistributable datasets.

OpenWRaw can represent mixed Waters triple-quadrupole acquisitions containing broad Q3 scans and MRM functions without requiring Waters SDK libraries.

## Data model

The TQ path keeps the two acquisition families distinct:

- **Broad Q3 scans** are decoded as dense m/z-intensity spectra from the low-resolution 6-byte packed representation.
- **MRM functions** are decoded from 4-byte packed intensities and emitted as one selected-reaction-monitoring chromatogram per Q1 -> Q3 transition.

By default, the mixed mzML writer therefore produces:

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

The canonical MRM chromatograms are unaffected by this choice.

## Optional MRM pseudo-MS2 compatibility stream

Canonical MRM data is chromatographic, not a native spectrum stream. OpenWRaw therefore keeps SRM chromatograms as the default and authoritative representation.

For spectrum-oriented downstream software, the mixed converter can **add** an explicit compatibility projection with `--mrm-pseudo-ms2`:

- transitions are grouped by Q1 precursor mass for each acquisition cycle;
- each Q1 group becomes one sparse MS2 spectrum;
- Q3 product masses form the m/z array;
- transition signals form the intensity array;
- Q1 is written as the precursor/selected m/z;
- the spectrum filter metadata explicitly identifies the record as an OpenWRaw pseudo-MS2 projection.

This option is additive. The original SRM chromatograms remain in `chromatogramList`, so the compatibility projection does not replace or destroy the canonical targeted representation.

When enabled, the mixed structure becomes:

```text
mzML
├── spectrumList
│   ├── broad Q3 scans
│   └── optional MRM pseudo-MS2 spectra
└── chromatogramList
    └── canonical MRM / SRM transition chromatograms
```

## Headless mixed conversion

From a source checkout, a format-faithful conversion is:

```bash
cargo run -p openwraw --example tq_mixed_to_mzml --release -- \
  path/to/bundle.raw output.mzML native-ms2
```

For an untargeted workflow that needs the broad Q3 scans presented as MS1:

```bash
cargo run -p openwraw --example tq_mixed_to_mzml --release -- \
  path/to/bundle.raw output.mzML pseudo-ms1
```

For a spectrum-oriented workflow that also needs MRM transitions surfaced as sparse MS2 records:

```bash
cargo run -p openwraw --example tq_mixed_to_mzml --release -- \
  path/to/bundle.raw output.mzML pseudo-ms1 --mrm-pseudo-ms2
```

Add `--indexed` to any of these commands for indexed mzML:

```bash
cargo run -p openwraw --example tq_mixed_to_mzml --release -- \
  path/to/bundle.raw output.mzML pseudo-ms1 --mrm-pseudo-ms2 --indexed
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

For canonical mixed mzML output:

```rust
use openwraw::tq_mixed_mzml::write_tq_mixed_mzml;
use openwraw::tq_mzml::TqQ3MzmlMode;

write_tq_mixed_mzml(
    "sample.raw",
    &mut output,
    TqQ3MzmlMode::PseudoMs1,
)?;
```

For the additive MRM pseudo-MS2 compatibility projection:

```rust
use openwraw::tq_mixed_mzml::write_tq_mixed_mzml_with_options;
use openwraw::tq_mrm_spectra::TqMrmSpectrumMode;
use openwraw::tq_mzml::TqQ3MzmlMode;

write_tq_mixed_mzml_with_options(
    "sample.raw",
    &mut output,
    TqQ3MzmlMode::PseudoMs1,
    TqMrmSpectrumMode::PseudoMs2,
)?;
```

## MZmine note

Modern MZmine's mzML importer parses `chromatogramList` together with precursor/product structures, so the canonical SRM representation is preserved at import. MZmine's dedicated MRM-to-scans processing workflow may require its MRM service depending on the distribution/license in use.

The optional `--mrm-pseudo-ms2` projection exists for cases where a downstream workflow primarily consumes spectrum lists. It is not required for archival conversion and should not be treated as a replacement for the canonical SRM chromatograms.

## Confidentiality and regression fixtures

The TQ implementation was developed using format-level observations that can be expressed independently of any one acquisition. Closed or confidential RAW data must **not** be committed as repository fixtures, copied into documentation, or reproduced as run-specific hexadecimal examples.

Repository tests for this path must use either:

1. synthetic binary records constructed by the tests themselves, or
2. independently redistributable public datasets whose licensing and provenance allow their use.

Private data may be used locally to verify that a generic decoder behaves correctly, but tests and documentation committed to the repository must not contain sample identifiers, private method details, real transition lists, chromatographic results, internal paths, serial numbers, or byte sequences copied from a closed run.

## Current limitations

The mixed converter currently focuses on the information required for structurally correct Q3 spectra and MRM chromatograms. Additional Waters side-file metadata such as compound labels and per-transition method parameters can be exposed later once their API and mzML representation are defined cleanly.

The legacy QTOF/IMS `Reader::open` path is intentionally unchanged. TQ support currently uses the dedicated TQ readers and conversion functions documented above.
