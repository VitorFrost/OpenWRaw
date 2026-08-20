//! Low-level helpers for Waters triple-quadrupole / low-resolution MassLynx RAW data.
//!
//! This module is deliberately format-focused. It does not change the existing
//! QTOF/IMS reader path and does not embed any private acquisition fixtures.
//! Tests use synthetic byte records only.

use crate::raw::data::Spectrum;
use crate::raw::header::FunctionCal;

/// Size of one `_FUNCTNS.INF` function record.
pub const FUNCTION_RECORD_SIZE: usize = 416;
/// Size of the legacy/simple `_FUNCnnn.IDX` record used by the TQ files
/// investigated for this format family.
pub const IDX22_STRIDE: usize = 22;
/// Number of Q1/Q3 slots available in a function descriptor.
pub const TRANSITION_SLOTS: usize = 32;
/// Mask for the pair-count field in the packed IDX word at +0x04.
///
/// Waters low-resolution readers use the lower 22 bits, not merely the lower
/// 16 bits. Keeping all 22 bits avoids truncating large scans.
pub const IDX_PAIR_COUNT_MASK: u32 = 0x003f_ffff;

/// Acquisition family recognized from the first byte of a TQ
/// `_FUNCTNS.INF` record.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TqFunctionKind {
    /// Multiple-reaction-monitoring function.
    Mrm,
    /// Broad MS2/Q3 scanning function.
    Q3Scan,
    /// Function code not currently classified by this module.
    Unknown(u8),
}

/// Per-function polarity encoded by the observed TQ function-type codes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TqPolarity {
    Positive,
    Negative,
}

/// Parsed static descriptor for one TQ function.
#[derive(Debug, Clone)]
pub struct TqFunctionDescriptor {
    pub function_type: u8,
    pub kind: TqFunctionKind,
    pub polarity: Option<TqPolarity>,
    /// Function-level set mass / precursor-context field at +0x018.
    pub set_mass_da: f32,
    /// Scan duration at +0x020.
    pub scan_time_s: f32,
    /// Lower acquisition bound at +0x0A0 for scan functions.
    pub mz_low: f32,
    /// Upper acquisition bound at +0x120 for scan functions.
    pub mz_high: f32,
    /// Per-slot dwell/timing values beginning at +0x020. Meaning is
    /// acquisition-mode dependent; primarily useful for MRM functions.
    pub dwell_or_timing: [f32; TRANSITION_SLOTS],
    /// Q1/parent mass slots beginning at +0x0A0 for MRM functions.
    pub q1_masses: [f32; TRANSITION_SLOTS],
    /// Q3/product mass slots beginning at +0x120 for MRM functions.
    pub q3_masses: [f32; TRANSITION_SLOTS],
}

impl TqFunctionDescriptor {
    /// Parse one exact 416-byte `_FUNCTNS.INF` record.
    pub fn from_record(record: &[u8]) -> crate::Result<Self> {
        if record.len() != FUNCTION_RECORD_SIZE {
            return Err(crate::Error::Parse(format!(
                "TQ function descriptor: expected {FUNCTION_RECORD_SIZE} bytes, got {}",
                record.len()
            )));
        }

        let function_type = record[0];
        let (kind, polarity) = classify_function_type(function_type);

        Ok(Self {
            function_type,
            kind,
            polarity,
            set_mass_da: crate::bytes::read_f32_le(record, 0x018)?,
            scan_time_s: crate::bytes::read_f32_le(record, 0x020)?,
            mz_low: crate::bytes::read_f32_le(record, 0x0a0)?,
            mz_high: crate::bytes::read_f32_le(record, 0x120)?,
            dwell_or_timing: read_f32_array32(record, 0x020)?,
            q1_masses: read_f32_array32(record, 0x0a0)?,
            q3_masses: read_f32_array32(record, 0x120)?,
        })
    }

    pub fn is_mrm(&self) -> bool {
        self.kind == TqFunctionKind::Mrm
    }

    pub fn is_q3_scan(&self) -> bool {
        self.kind == TqFunctionKind::Q3Scan
    }
}

/// Classify TQ function-type bytes observed across positive and negative
/// acquisition modes.
///
/// Exact matches are used intentionally rather than treating bit 0x20 as a
/// universal polarity bit for every Waters instrument family.
pub fn classify_function_type(code: u8) -> (TqFunctionKind, Option<TqPolarity>) {
    match code {
        0x09 => (TqFunctionKind::Mrm, Some(TqPolarity::Positive)),
        0x29 => (TqFunctionKind::Mrm, Some(TqPolarity::Negative)),
        0x0b => (TqFunctionKind::Q3Scan, Some(TqPolarity::Positive)),
        0x2b => (TqFunctionKind::Q3Scan, Some(TqPolarity::Negative)),
        other => (TqFunctionKind::Unknown(other), None),
    }
}

fn read_f32_array32(data: &[u8], offset: usize) -> crate::Result<[f32; TRANSITION_SLOTS]> {
    let mut out = [0.0_f32; TRANSITION_SLOTS];
    for (i, value) in out.iter_mut().enumerate() {
        *value = crate::bytes::read_f32_le(data, offset + i * 4)?;
    }
    Ok(out)
}

/// One 22-byte `_FUNCnnn.IDX` record for low-resolution/TQ data.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TqIndexRecord {
    pub dat_offset: u32,
    pub packed: u32,
    pub pair_count: u32,
    pub retention_time_min: f32,
}

/// Parse a complete 22-byte-stride IDX file.
pub fn parse_idx22(data: &[u8]) -> crate::Result<Vec<TqIndexRecord>> {
    if data.len() % IDX22_STRIDE != 0 {
        return Err(crate::Error::Parse(format!(
            "TQ IDX: size {} is not a multiple of {IDX22_STRIDE}",
            data.len()
        )));
    }

    let mut out = Vec::with_capacity(data.len() / IDX22_STRIDE);
    for record in data.chunks_exact(IDX22_STRIDE) {
        let dat_offset = crate::bytes::read_u32_le(record, 0x00)?;
        let packed = crate::bytes::read_u32_le(record, 0x04)?;
        let retention_time_min = crate::bytes::read_f32_le(record, 0x0c)?;
        out.push(TqIndexRecord {
            dat_offset,
            packed,
            pair_count: packed & IDX_PAIR_COUNT_MASK,
            retention_time_min,
        });
    }
    Ok(out)
}

/// Infer the DAT record width from the final non-empty IDX entry.
///
/// This mirrors the format-intrinsic relationship
/// `(dat_size - final_offset) / final_pair_count` and avoids assuming that a
/// 22-byte IDX always implies one particular DAT encoding.
pub fn infer_bytes_per_pair(dat_size: u64, index: &[TqIndexRecord]) -> crate::Result<usize> {
    let last = index
        .iter()
        .rev()
        .find(|r| r.pair_count != 0)
        .ok_or_else(|| crate::Error::Parse("TQ IDX: no non-empty records".to_owned()))?;

    let offset = last.dat_offset as u64;
    if offset > dat_size {
        return Err(crate::Error::Parse(format!(
            "TQ IDX: final DAT offset {offset} exceeds DAT size {dat_size}"
        )));
    }

    let remaining = dat_size - offset;
    let count = last.pair_count as u64;
    if remaining % count != 0 {
        return Err(crate::Error::Parse(format!(
            "TQ DAT: trailing byte count {remaining} is not divisible by pair count {count}"
        )));
    }

    let width = (remaining / count) as usize;
    if !matches!(width, 2 | 4 | 6 | 8) {
        return Err(crate::Error::Parse(format!(
            "TQ DAT: unsupported inferred record width {width} bytes"
        )));
    }
    Ok(width)
}

/// Return the byte range for one scan/cycle using a known record width.
pub fn scan_byte_range(
    index: &[TqIndexRecord],
    scan_index: usize,
    bytes_per_pair: usize,
    dat_size: u64,
) -> crate::Result<std::ops::Range<usize>> {
    let rec = index
        .get(scan_index)
        .ok_or_else(|| crate::Error::Parse(format!("TQ scan index {scan_index} out of range")))?;
    let start = rec.dat_offset as u64;
    let length = (rec.pair_count as u64)
        .checked_mul(bytes_per_pair as u64)
        .ok_or_else(|| crate::Error::Parse("TQ scan byte length overflow".to_owned()))?;
    let end = start
        .checked_add(length)
        .ok_or_else(|| crate::Error::Parse("TQ scan byte range overflow".to_owned()))?;
    if end > dat_size {
        return Err(crate::Error::Parse(format!(
            "TQ scan byte range {start}..{end} exceeds DAT size {dat_size}"
        )));
    }
    Ok(start as usize..end as usize)
}

/// Decode the direct 6-byte low-resolution m/z-intensity representation.
///
/// Record layout:
/// - bytes 0..2: signed 16-bit intensity base
/// - bytes 2..6: packed mass base/exponent plus intensity exponent
///
/// Unlike OpenWRaw's TOF Encoding A, this representation stores m/z directly;
/// it has no per-scan TOF sentinel and needs no TOF flight-path geometry.
pub fn decode_direct6(
    scan_bytes: &[u8],
    calibration: Option<&FunctionCal>,
) -> crate::Result<Spectrum> {
    if scan_bytes.len() % 6 != 0 {
        return Err(crate::Error::Parse(format!(
            "TQ direct-6: scan size {} is not a multiple of 6",
            scan_bytes.len()
        )));
    }

    let n = scan_bytes.len() / 6;
    let mut spectrum = Spectrum {
        mz: Vec::with_capacity(n),
        intensity: Vec::with_capacity(n),
    };

    for record in scan_bytes.chunks_exact(6) {
        let intensity_base = i16::from_le_bytes([record[0], record[1]]) as f64;
        let packed = u32::from_le_bytes([record[2], record[3], record[4], record[5]]);

        let mass_base = packed >> 9;
        let mass_power = ((packed & 0x01f0) >> 4) as i32 - 23;
        let raw_mz = (mass_base as f64) * 2.0_f64.powi(mass_power);
        let mz = calibration
            .filter(|cal| !cal.coeffs.is_empty())
            .map(|cal| apply_direct_mass_calibration(raw_mz, cal))
            .unwrap_or(raw_mz);

        let intensity_power = (packed & 0x0f) as i32;
        let intensity = intensity_base * 4.0_f64.powi(intensity_power);

        spectrum.mz.push(mz);
        spectrum.intensity.push(intensity as f32);
    }

    Ok(spectrum)
}

/// Apply the calibration coefficients as a polynomial directly in the mass
/// domain. Low-resolution packed-mass data uses this convention independently
/// of the TOF-oriented T0/T1 interpretation used by the QTOF decoder path.
fn apply_direct_mass_calibration(raw_mz: f64, calibration: &FunctionCal) -> f64 {
    if calibration.coeffs.is_empty() {
        return raw_mz;
    }
    let mut result = 0.0_f64;
    for &coefficient in calibration.coeffs.iter().rev() {
        result = result * raw_mz + coefficient;
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::raw::header::{CalType, FunctionCal};

    fn direct6_record(
        mass_base: u32,
        mass_power_field: u32,
        intensity_base: i16,
        intensity_power: u32,
    ) -> [u8; 6] {
        let packed = (mass_base << 9) | ((mass_power_field & 0x1f) << 4) | (intensity_power & 0x0f);
        let mut record = [0_u8; 6];
        record[0..2].copy_from_slice(&intensity_base.to_le_bytes());
        record[2..6].copy_from_slice(&packed.to_le_bytes());
        record
    }

    #[test]
    fn classifies_tq_function_codes() {
        assert_eq!(
            classify_function_type(0x09),
            (TqFunctionKind::Mrm, Some(TqPolarity::Positive))
        );
        assert_eq!(
            classify_function_type(0x29),
            (TqFunctionKind::Mrm, Some(TqPolarity::Negative))
        );
        assert_eq!(
            classify_function_type(0x0b),
            (TqFunctionKind::Q3Scan, Some(TqPolarity::Positive))
        );
        assert_eq!(
            classify_function_type(0x2b),
            (TqFunctionKind::Q3Scan, Some(TqPolarity::Negative))
        );
        assert_eq!(
            classify_function_type(0x55),
            (TqFunctionKind::Unknown(0x55), None)
        );
    }

    #[test]
    fn parses_synthetic_q3_function_descriptor() {
        let mut record = [0_u8; FUNCTION_RECORD_SIZE];
        record[0] = 0x0b;
        record[0x018..0x01c].copy_from_slice(&40.0_f32.to_le_bytes());
        record[0x020..0x024].copy_from_slice(&0.5_f32.to_le_bytes());
        record[0x0a0..0x0a4].copy_from_slice(&75.0_f32.to_le_bytes());
        record[0x120..0x124].copy_from_slice(&900.0_f32.to_le_bytes());

        let descriptor = TqFunctionDescriptor::from_record(&record).unwrap();
        assert!(descriptor.is_q3_scan());
        assert_eq!(descriptor.polarity, Some(TqPolarity::Positive));
        assert!((descriptor.set_mass_da - 40.0).abs() < f32::EPSILON);
        assert!((descriptor.mz_low - 75.0).abs() < f32::EPSILON);
        assert!((descriptor.mz_high - 900.0).abs() < f32::EPSILON);
    }

    #[test]
    fn idx_parser_keeps_22_bit_pair_count() {
        let pair_count = 0x02_aaaa_u32;
        let packed = 0x1800_0000_u32 | pair_count;
        let mut record = [0_u8; IDX22_STRIDE];
        record[0..4].copy_from_slice(&1234_u32.to_le_bytes());
        record[4..8].copy_from_slice(&packed.to_le_bytes());
        record[12..16].copy_from_slice(&1.25_f32.to_le_bytes());

        let parsed = parse_idx22(&record).unwrap();
        assert_eq!(parsed.len(), 1);
        assert_eq!(parsed[0].dat_offset, 1234);
        assert_eq!(parsed[0].pair_count, pair_count);
        assert!((parsed[0].retention_time_min - 1.25).abs() < f32::EPSILON);
    }

    #[test]
    fn infers_record_width_from_last_nonempty_scan() {
        let index = vec![
            TqIndexRecord {
                dat_offset: 0,
                packed: 2,
                pair_count: 2,
                retention_time_min: 0.0,
            },
            TqIndexRecord {
                dat_offset: 12,
                packed: 3,
                pair_count: 3,
                retention_time_min: 0.1,
            },
        ];
        assert_eq!(infer_bytes_per_pair(30, &index).unwrap(), 6);
    }

    #[test]
    fn direct6_decodes_mass_and_intensity() {
        // mass_power_field=23 makes the mass scale 2^(23-23)=1.
        let record = direct6_record(100, 23, 42, 0);
        let spectrum = decode_direct6(&record, None).unwrap();
        assert_eq!(spectrum.mz, vec![100.0]);
        assert_eq!(spectrum.intensity, vec![42.0]);
    }

    #[test]
    fn direct6_decodes_fractional_mass_and_scaled_intensity() {
        // mass_power_field=22 -> scale 0.5; 201 * 0.5 = 100.5.
        // intensity exponent 2 -> base * 4^2 = base * 16.
        let record = direct6_record(201, 22, 10, 2);
        let spectrum = decode_direct6(&record, None).unwrap();
        assert!((spectrum.mz[0] - 100.5).abs() < 1e-12);
        assert_eq!(spectrum.intensity[0], 160.0);
    }

    #[test]
    fn direct6_applies_mass_domain_polynomial_even_for_t0_label() {
        let record = direct6_record(100, 23, 1, 0);
        let calibration = FunctionCal {
            coeffs: vec![1.0, 2.0],
            cal_type: CalType::T0,
        };
        let spectrum = decode_direct6(&record, Some(&calibration)).unwrap();
        assert!((spectrum.mz[0] - 201.0).abs() < 1e-12);
    }

    #[test]
    fn scan_range_is_bounds_checked() {
        let index = vec![TqIndexRecord {
            dat_offset: 12,
            packed: 3,
            pair_count: 3,
            retention_time_min: 0.1,
        }];
        assert_eq!(scan_byte_range(&index, 0, 6, 30).unwrap(), 12..30);
        assert!(scan_byte_range(&index, 0, 6, 29).is_err());
    }
}
