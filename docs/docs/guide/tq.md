# Triple-quadrupole / low-resolution Waters RAW

> Experimental. This path is intentionally separate from the established QTOF/IMS `Reader` while TQ support is validated across independently redistributable datasets.

OpenWRaw can represent mixed Waters triple-quadrupole acquisitions containing broad Q3 scans and MRM functions without requiring Waters SDK libraries.

## Validation status

The broad-Q3 direct-6 decoder has been exercised against a complete non-public TQ Q3 function in addition to the repository's synthetic tests. Without publishing any private fixture data, the validation confirmed exact IDX-to-DAT coverage, successful decoding of every scan, monotonic m/z ordering, and agreement between decoded intensity sums and the IDX `+0x08` scan-signal statistic at relative error on the order of `10^-6`.

The direct packed m/z equation was also compared against Rainbow's independent public Waters implementation. The uncalibrated m/z values agreed point-for-point; after calibration, the remaining difference was below `2.3e-4` Da (`0.22` ppm) and is explained by Rainbow's `float32` arithmetic versus OpenWRaw's `f64` polynomial evaluation. Vendor/ProteoWizard point-for-point confirmation of the calibrated m/z axis is still pending.

The MRM DAT4 path has now been structurally checked across the private acquisition's positive- and negative-polarity MRM functions and different transition counts. IDX cycle geometry and DAT length close exactly. A repeatable encoding-specific relationship was also found: decoded DAT4 intensity sums are approximately twice the IDX `+0x08` statistic. The same 2:1 relationship occurs in Rainbow's public Waters TQ fixture, so OpenWRaw deliberately does **not** introduce a factor-of-two correction merely to make DAT4 agree with that IDX field. Absolute DAT4 intensity scale still requires vendor-reference confirmation.

The complete private mixed TQ bundle has also been exercised through the mixed source and PSI mzML writers in native-MS2, pseudo-MS1, and additive MRM pseudo-MS2 modes. The broad-Q3 binary payload is unchanged when the MRM compatibility projection is enabled, acquisition-time ordering remains monotonic, canonical SRM chromatograms remain present, and both plain and indexed outputs pass the bundled PSI mzML schemas. Indexed output has additionally been checked for spectrum/chromatogram offsets and file-checksum consistency. These validation statements intentionally omit run-specific counts and acquisition details.

TQ support therefore remains **experimental**, especially outside the observed instrument/function families.

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

Broad Q3 functions are merged by retention time rather than emitted as separate function-sized blocks. This preserves acquisition chronology in polarity-switching or otherwise interleaved runs while each spectrum retains its Waters `function/process/scan` identity.

## Q3 semantics

Waters triple-quadrupole broad scans may be recorded as MS2/Q3 scans even when the resulting broad mass spectrum is intended for untargeted downstream processing. OpenWRaw therefore requires an explicit choice rather than silently changing the MS level.

### `native-ms2`

Preserves the vendor acquisition semantics. Broad Q3 scans remain MS2 and retain valid non-zero function-level precursor/set-mass context when present.

Use this for format-faithful interchange.

### `pseudo-ms1`

Projects the broad Q3 spectrum to MS1 for downstream tools that require an MS1 survey stream. The transformation is explicit and is written as an `openwraw.projection` `userParam`; OpenWRaw does **not** misuse the Thermo-specific PSI `filter string` term for this annotation.

The canonical MRM chromatograms are unaffected by this choice.

## Optional MRM pseudo-MS2 compatibility stream

Canonical MRM data is chromatographic, not a native spectrum stream. OpenWRaw therefore keeps SRM chromatograms as the default and authoritative representation.

For spectrum-oriented downstream software, the mixed converter can **add** an explicit compatibility projection with `--mrm-pseudo-ms2`:

- transitions are grouped by Q1 precursor mass for each acquisition cycle;
- each Q1 group becomes one sparse MS2 spectrum;
- Q3 product masses form the m/z array and are sorted while retaining their paired intensities;
- transition signals form the intensity array;
- Q1 is written as the precursor/selected m/z;
- `openwraw.projection=pseudo-ms2-from-mrm` records that the spectrum is generated rather than native.

This option is additive. The original SRM chromatograms remain in `chromatogramList`, so the compatibility projection does not replace or destroy the canonical targeted representation.

When pseudo-MS2 compatibility records are enabled, they are merged with broad Q3 spectra by retention time and the final spectrum indices are reassigned contiguously after that merge.

When enabled, the mixed structure becomes:

```text
mzML
├── spectrumList
│   ├── broad Q3 scans
│   └── optional MRM pseudo-MS2 spectra
└── chromatogramList
    └── canonical MRM / SRM transition chromatograms
```

## PSI mzML semantics

The TQ conversion path applies an additional PSI-oriented serialization pass after the generic OpenMassSpecCore 1.5.0 writer. This is deliberately scoped to TQ output so the established QTOF/IMS writer behavior is not changed by this experimental feature.

The correction layer currently ensures that:

- `fileContent` advertises only MS1/MSn spectrum types actually emitted;
- Waters source files are represented individually with per-file SHA-1 provenance while the actual private source names remain separate from generated mzML identifiers;
- an unknown MS2/SRM dissociation mechanism is represented by the generic PSI `MS:1000044` `dissociation method` term rather than an empty activation or invented CID assertion;
- projection provenance is carried in `userParam` values instead of `MS:1000512 filter string`;
- intensity arrays and spectrum-level mass/intensity summary values carry appropriate PSI units;
- the instrument configuration contains a PSI `componentList` with a conservative ionization source, Q1 quadrupole, Q3 quadrupole, and detector. Generic parent terms are used for source/detector type until more specific native metadata is decoded;
- Q3 `lowest/highest observed m/z` remain based on the decoded spectrum, while `scanWindow` uses the programmed lower/upper acquisition bounds from `_FUNCTNS.INF`;
- the mzML `run` identifier is normalized to a conservative XML `xs:ID`-safe value when a RAW directory name starts with a number or contains unsupported punctuation, without changing the reconstructed source-file provenance;
- indexed mzML offsets are calculated only after these semantic corrections, and `fileChecksum` is recomputed over the corrected indexed document according to the indexed-mzML checksum convention.

The TQ serializer is tied to the exact `openmassspec-core` version in `Cargo.lock`; its string-level corrections must be reviewed whenever that dependency's mzML writer changes.

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

**Generated mzML from a confidential RAW is itself confidential.** It contains source-derived spectral/chromatographic data plus per-file source provenance and cryptographic checksums. Do not publish a converted private mzML merely because the OpenWRaw source code and synthetic tests are safe to publish.

## Current limitations

The mixed converter currently focuses on the information required for structurally and semantically defensible Q3 spectra and MRM chromatograms. Additional Waters side-file metadata such as compound labels and per-transition method parameters can be exposed later once their API and mzML representation are defined cleanly.

`openmassspec-core` 1.5 exposes an infallible `SpectrumSource` iterator. OpenWRaw validates static Q3 IDX/DAT layout before iteration and deliberately fails fast if a Q3 DAT read later becomes impossible (for example because the source files are removed or modified during conversion), rather than silently omitting a spectrum. A future fallible iterator in the core API would allow this environmental failure to be returned as a normal `Result` instead of a hard failure.

MRM functions are currently required to have a constant non-zero transition count across their populated acquisition cycles. Zero-point cycles are skipped. Scheduled/dynamic MRM functions that change the active transition count within one function are rejected rather than guessed, because the current format model does not yet establish an unambiguous per-cycle mapping from the changing DAT channel set back to the 32 Q1/Q3 descriptor slots.

The current MRM transition model also assumes the active Q1/Q3 descriptor entries occupy the leading transition slots. This is validated for the format family used to develop the reader but has not yet been generalized to sparse/non-contiguous descriptor-slot schedules.

DAT4 absolute signal scaling has not yet been compared against a Waters SDK / MassLynx / ProteoWizard reference. The observed 2:1 relationship between decoded DAT4 sums and IDX `+0x08` is reproducible across both the private TQ dataset and a public Rainbow TQ fixture, but `+0x08` itself does not have a universal TIC scale across MassLynx encodings and is not used to rescale MRM output.

Specific ion-source and detector component types are not yet decoded from the TQ RAW path, so the mzML component list uses PSI parent terms for those two components while representing Q1 and Q3 explicitly as quadrupoles.

The legacy QTOF/IMS `Reader::open` path is intentionally unchanged. TQ support currently uses the dedicated TQ readers and conversion functions documented above.
