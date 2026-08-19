# Triple-quadrupole / low-resolution MassLynx variants

## Status: Partially decoded, privately validated

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

OpenWRaw currently masks only the low 16 bits for Variant A. Expanding this mask is recommended before supporting additional low-resolution files.

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

For canonical mzML, MRM data should preferably be exported as SRM chromatograms rather than fabricated full spectra. A compatibility option analogous to `srmAsSpectra` may be added separately.

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

This should be implemented as a separate encoding variant, for example `LowResPacked6`, rather than modifying the existing TOF `Encoding A`.

## `_extern.inf` implications

Current OpenWRaw requires TOF geometry fields such as `Lteff`, `Veff`, and a pusher interval during `Reader::open()`.

Those fields are not intrinsic requirements of MassLynx RAW and are not required to decode the direct low-resolution TQ format described above.

Recommended model:

```text
ExternInf
└── optional TOF geometry
    ├── Lteff
    ├── Veff
    └── pusher interval
```

TOF decoders should require the geometry only when the selected DAT encoding actually needs it.

### Polarity is per function

Mixed-polarity triple-quadrupole methods may contain both positive and negative functions in the same RAW bundle. A single run-level polarity is therefore insufficient.

Recommended change:

```text
ExternFunction
└── polarity: Option<Polarity>
```

The mzML exporter should resolve polarity using the current function index.

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

## Recommended implementation changes

1. Add TQ function-kind decoding based on `function_type`.
2. Treat TOF geometry in `_extern.inf` as optional.
3. Expand the Variant-A/22-byte pair-count mask to 22 bits.
4. Detect DAT width using function type plus DAT/IDX consistency, not IDX stride alone.
5. Add a `LowResPacked6` decoder for direct m/z/intensity records.
6. Add a 4-byte MRM intensity decoder.
7. Parse Q1 and Q3 transition arrays from `_FUNCTNS.INF`.
8. Store polarity per function.
9. Export MRM as chromatograms and broad Q3 scans as native MS2 by default.
10. Add an explicit pseudo-MS1 compatibility mode only when requested.

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
- full `.EE` and `.CMP` generalization across instrument generations.

The parser should therefore preserve unknown raw fields where practical and avoid over-generalizing from a single instrument family.
