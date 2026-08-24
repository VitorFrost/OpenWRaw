# Triple-quadrupole / low-resolution MassLynx variants

## Status: Q3 direct-6 validated on private TQ data; MRM DAT4 structurally cross-validated; broader TQ support remains experimental

This note documents format-level behavior observed in a non-public Waters triple-quadrupole MassLynx RAW dataset and cross-checked against public implementations where possible.

**Confidentiality rule:** no sample identifiers, acquisition dates, instrument serial numbers, internal file-system paths, method names, compound names, transition lists, chromatographic results, raw byte excerpts, or other run-specific values from the non-public dataset are included here. The private dataset must not be committed to this repository or distributed as a test fixture.

The goal is to capture only reusable binary-format knowledge needed to extend OpenWRaw safely.

## Scope

The validated RAW bundle uses the classic MassLynx directory layout:

- `_HEADER.TXT`
- `_FUNCTNS.INF`
- `_extern.inf`
- `_FUNCnnn.IDX`
- `_FUNCnnn.DAT`
- optional per-function files such as `_FUNCnnn.STS`, `.EE`, and `.CMP`

The acquisition contains both targeted MRM functions and broad MS2/Q3 scan functions. This is a useful compatibility case because the current OpenWRaw corpus is dominated by TOF/QTOF and IMS-TOF instruments.

### Private Q3 validation result

The 6-byte broad-Q3 path has now been exercised end-to-end against a complete non-public TQ Q3 function, without committing any source-derived fixture data. The validation used only format invariants and aggregate numerical checks:

- every IDX scan range resolved inside the paired DAT file;
- the scan ranges covered the DAT payload exactly, with no gaps, overlaps, or trailing bytes;
- every scan decoded successfully through `TqReader`;
- decoded m/z values were monotonically ordered within every scan;
- the sum of decoded point intensities reproduced the scan-signal statistic stored at IDX `+0x08` to relative error on the order of `10^-6`;
- the per-function calibration polynomial from `_HEADER.TXT` was applied in the direct mass domain by the reader.

These checks strongly validate the 22-byte IDX interpretation, 22-bit pair-count mask, 6-byte direct packed representation, intensity decoder, and scan slicing logic for this TQ family. The `+0x08` relationship is an empirical invariant for the observed direct-6 Q3 data and must not be generalized as a universal TIC definition for 22-byte IDX records. Absolute calibrated m/z accuracy has **not yet** been independently compared point-for-point against a vendor-generated or ProteoWizard reference export, so that remains a separate validation step.

### Independent direct-6 implementation check

The packed-mass equation and the 22-bit IDX pair-count rule were independently compared against Rainbow's public Waters parser. Across the complete private Q3 function, the uncalibrated packed m/z values were identical point-for-point to the Rainbow equation. Applying the same calibration polynomial differed only by numerical precision: Rainbow accumulates its calibrated axis in `float32`, whereas OpenWRaw evaluates the polynomial in `f64`. The largest observed difference was below `2.3e-4` Da (below `0.22` ppm).

This is independent implementation agreement, not a substitute for a vendor-reference export.

## `_FUNCTNS.INF`

The common 416-byte-per-function structure remains valid for the triple-quadrupole dataset.

### Function type codes observed

The following values were observed and validated against the paired DAT/IDX behavior:

| `+0x000` | Interpretation | Polarity |
|---:|---|---|
| `0x09` | MRM | positive |
| `0x29` | MRM | negative |
| `0x0B` | MS2 / Q3 broad scan | positive |
| `0x2B` | MS2 / Q3 broad scan | negative |

The difference `0x20` is consistent with a polarity flag in this instrument family. This should currently be treated as an **observed TQ-family rule**, not as a universal MassLynx guarantee across all instruments.

`scan_subtype` at `+0x001` must not be interpreted without considering `function_type`: values already used in the TOF/QTOF corpus may occur with different function families.

### MRM-specific arrays

For MRM function records, three 32-slot `f32 LE` arrays were validated:

| Offset range | Layout | Meaning |
|---|---|---|
| `0x020..0x09F` | 32 × `f32` | transition timing/dwell-related values |
| `0x0A0..0x11F` | 32 × `f32` | Q1 / precursor masses |
| `0x120..0x19F` | 32 × `f32` | Q3 / product masses |

The Q3 array location at `+0x120` is independently consistent with Rainbow's public MassLynx parser.

For broad MS2/Q3 scan functions, the same record positions are interpreted differently:

- `+0x018`: function-level Set Mass / precursor-context value
- `+0x020`: scan time
- `+0x0A0`: scan lower m/z bound
- `+0x120`: scan upper m/z bound

The Set Mass value is separate from the scan window. It must not be treated as the lower scan mass.

## `_FUNCnnn.IDX`: 22-byte low-resolution variant

The triple-quadrupole functions use the existing 22-byte IDX stride.

Fields validated for this family:

| Offset | Type | Meaning |
|---:|---|---|
| `+0x00` | `u32 LE` | DAT byte offset for the acquisition cycle / scan |
| `+0x04` | `u32 LE` | packed format/count field |
| `+0x0C` | `f32 LE` | retention time in minutes |

The low 22 bits of the packed field are used as the pair/value count in public Waters parsers:

```text
pair_count = packed & 0x003F_FFFF
```

The dedicated TQ parser uses the full 22-bit count field. The established TOF Variant-A path is kept separate so its legacy behavior is not changed solely from observations in the TQ family.

### Do not infer DAT encoding from IDX stride alone

A 22-byte IDX does **not** imply a single 6-byte DAT encoding.

In the validated TQ family:

- MRM functions use 4-byte packed intensity values;
- broad low-resolution scan functions use 6-byte packed `(m/z, intensity)` records.

Therefore encoding detection should combine:

1. IDX stride,
2. function type,
3. DAT size,
4. DAT offset/count consistency.

A robust generic bytes-per-pair check is:

```text
bytes_per_pair = (dat_size - final_dat_offset) / final_pair_count
```

when the final non-empty IDX entry is well formed.

## 4-byte MRM DAT encoding

MRM DAT files contain intensity values only. Q1/Q3 coordinates are supplied by `_FUNCTNS.INF`.

For a function with `M` active transitions and `N` acquisition cycles:

```text
DAT size = N × M × 4 bytes
```

Each 4-byte little-endian value decodes as:

```text
raw   = u32_le
power = raw >> 22
base  = raw & 0x001F_FFFF

intensity = (base / 1024) × 2^(power - 10)
```

The output data model is therefore conceptually:

```text
cycle 1: I(transition 1), I(transition 2), ...
cycle 2: I(transition 1), I(transition 2), ...
...
```

The transition m/z values must not be reconstructed from DAT because they are not stored there.

OpenWRaw exports this data canonically as SRM chromatograms rather than fabricating native full spectra. An optional additive pseudo-MS2 projection is available for spectrum-oriented compatibility workflows while preserving the SRM chromatograms as the authoritative representation.

### MRM DAT4 validation status

The 4-byte path has been checked across all MRM functions in the private mixed-polarity TQ acquisition. The functions covered both polarities and different active transition counts. For every tested function:

- the IDX pair count was constant within the function and matched the active leading Q1/Q3 descriptor slots;
- `cycle_count × transition_count × 4` matched the DAT size exactly;
- every IDX-derived DAT range was in bounds, with no overlaps or trailing bytes;
- decoding with the public DAT4 equation produced finite non-negative signal values.

An additional invariant was observed at IDX `+0x08`: for every tested MRM cycle, the sum of decoded DAT4 transition intensities was approximately **two times** the `f32` value stored at `+0x08`, with relative residuals on the order of `10^-7`. The same 2:1 relationship was independently reproduced on Rainbow's public Waters TQ fixture, which uses the same 4-byte encoding.

This cross-dataset result is important because it shows that IDX `+0x08` is **not a canonical TIC field with one universal scale across low-resolution encodings**. In direct-6 Q3 data it matches the decoded intensity sum, while in the observed DAT4 MRM data it tracks half of that sum. OpenWRaw therefore does **not** rescale DAT4 intensities merely to force equality with IDX `+0x08`.

The absolute intensity scale of DAT4 should remain marked provisional until compared against a vendor SDK / MassLynx / ProteoWizard reference. The current decoder intentionally follows the independently published Rainbow equation rather than introducing an unsupported factor-of-two correction.

## 6-byte direct low-resolution scan encoding

Broad TQ scan functions also use 6-byte records, but they are **not** the same representation as OpenWRaw's current TOF `Encoding A`.

OpenWRaw `Encoding A` expects a sentinel-anchored TOF-bin representation whose m/z axis is reconstructed from TOF geometry. The TQ low-resolution representation stores an m/z/intensity pair directly in each six-byte record and does not require TOF flight-path parameters.

A compatible public decoder is present in Rainbow. The record can be interpreted using a signed 16-bit intensity base and a packed 32-bit mass/intensity descriptor:

```text
base_value = i16_le(bytes[0..2])
raw        = u32_le(bytes[2..6])

value_power = raw & 0x0F
mass_power  = (raw & 0x1F0) >> 4
mass_base   = raw >> 9

mz = mass_base × 2^(mass_power - 23)
intensity = base_value × 4^(value_power)
```

Per-function calibration from `_HEADER.TXT` may then be applied to the decoded m/z values.

The TQ path implements this as a separate direct-6 decoder (`decode_direct6`) rather than modifying the existing TOF `Encoding A`.

## `_extern.inf` implications

The established QTOF/IMS `Reader::open()` path requires TOF geometry fields such as `Lteff`, `Veff`, and a pusher interval. Those fields are not intrinsic requirements of MassLynx RAW and are not required to decode the direct low-resolution TQ format described above.

The experimental TQ implementation therefore uses a dedicated `TqReader` that does not require `_extern.inf` TOF geometry. This avoids weakening or changing the validated TOF reader while TQ support is still being generalized. A future shared reader model may make analyzer-specific geometry optional and require it only in decoders that actually use TOF coordinates.

### Polarity is per function

Mixed-polarity triple-quadrupole methods may contain both positive and negative functions in the same RAW bundle. A single run-level polarity is therefore insufficient.

The TQ function descriptor therefore stores polarity per function, and the TQ mzML exporter resolves polarity from the current function rather than from a single run-level value. This behavior is scoped to the observed TQ function codes and is not assumed to be universal across Waters families.

## mzML semantics

The native acquisition semantics should be preserved by default.

### MRM

Export as SRM/MRM chromatograms with Q1/Q3 transition metadata.

### Broad Q3/MS2 scan

Export as an MS2 spectrum with:

- the true scan window from the function descriptor;
- the function-level Set Mass kept as separate precursor/context metadata when supported;
- decoded direct m/z/intensity pairs from the low-resolution 6-byte DAT representation.

For downstream tools that require MS1-like survey data, an **explicit opt-in transformation** may relabel broad Q3/MS2 scans as pseudo-MS1. This must be recorded as a computational reinterpretation, not presented as native MS1 acquisition.

## Implementation status

The experimental TQ branch now implements the core items originally identified during format analysis:

- TQ function-kind decoding from observed `function_type` values;
- a dedicated `TqReader` path that does not require TOF flight-path geometry;
- the 22-bit IDX pair-count mask;
- DAT record-width inference from IDX/DAT consistency;
- a 6-byte direct low-resolution m/z/intensity decoder for broad Q3 scans;
- a 4-byte packed-intensity decoder for MRM cycles;
- Q1/Q3 transition arrays from `_FUNCTNS.INF`;
- per-function TQ polarity;
- native MS2 export for broad Q3 scans;
- canonical SRM chromatogram export for MRM;
- explicit pseudo-MS1 projection for broad Q3 scans;
- optional additive pseudo-MS2 projection for MRM compatibility workflows.

The established QTOF/IMS `Reader` remains unchanged while this TQ path is validated across additional independent datasets. Generalizing the common reader architecture can be considered after the format families are better covered by corpus data.

## Testing and confidentiality policy

The non-public RAW dataset used to establish these observations must remain outside the repository.

Tests added to OpenWRaw should use only:

- synthetic 416-byte function descriptors;
- synthetic 22-byte IDX records;
- synthetic 4-byte intensity values;
- synthetic 6-byte direct m/z/intensity records;
- public datasets that are independently redistributable.

Do **not** copy into tests, comments, fixtures, documentation, issues, pull requests, CI artifacts, or commit messages any sample-specific metadata or byte sequences taken from a closed dataset.

When a regression test must reproduce a structural edge case found in private data, construct a minimal artificial record that exercises the same parser behavior without preserving the original numeric values or acquisition content.

## Open questions

The following should remain marked as provisional until confirmed across additional independent TQ-family files:

- whether the `0x20` polarity relationship applies to all Waters triple-quadrupole generations;
- the complete enumeration of TQ `function_type` values;
- exact semantics of every non-zero `_FUNCTNS.INF` field for MRM and broad scan functions;
- whether the function-level Set Mass in broad Q3 scans always corresponds to a meaningful Q1 setting or may sometimes be acquisition-software bookkeeping;
- independent point-for-point confirmation of the calibrated Q3 m/z axis against a vendor-generated or ProteoWizard reference export;
- absolute DAT4 MRM intensity-scale confirmation against a vendor-generated or ProteoWizard reference export;
- full `.EE` and `.CMP` generalization across instrument generations.

The parser should therefore preserve unknown raw fields where practical and avoid over-generalizing from a single instrument family.
