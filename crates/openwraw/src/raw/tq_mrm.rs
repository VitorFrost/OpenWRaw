//! Waters triple-quadrupole MRM decoding for MassLynx RAW bundles.
//!
//! MRM functions in this format family use a 22-byte `_FUNCnnn.IDX` stream
//! and a 4-byte packed intensity value per transition per acquisition cycle.
//! Q1/Q3 transition masses are read from the corresponding 416-byte
//! `_FUNCTNS.INF` record; the DAT payload therefore only needs to store signal.
//!
//! No compound-name or method-specific side-file metadata is required here.

use std::fs;
use std::path::{Path, PathBuf};

use crate::raw::tq::{
    infer_bytes_per_pair, parse_idx22, scan_byte_range, TqFunctionDescriptor, TqFunctionKind,
    TqIndexRecord, TqPolarity, FUNCTION_RECORD_SIZE, TRANSITION_SLOTS,
};

/// One monitored Q1 -> Q3 transition.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TqMrmTransition {
    pub slot: usize,
    pub precursor_mz: f64,
    pub product_mz: f64,
}

/// One MRM function ready for decoding.
#[derive(Debug, Clone)]
pub struct TqMrmFunctionEntry {
    pub index: u32,
    pub descriptor: TqFunctionDescriptor,
    pub scan_index: Vec<TqIndexRecord>,
    pub dat_path: PathBuf,
    pub dat_size: u64,
    pub bytes_per_pair: usize,
    pub transitions: Vec<TqMrmTransition>,
}

impl TqMrmFunctionEntry {
    pub fn cycle_count(&self) -> usize {
        self.scan_index.len()
    }
}

/// Canonical chromatographic representation of one monitored reaction.
#[derive(Debug, Clone)]
pub struct TqMrmChromatogram {
    pub function_index: u32,
    pub transition_index: usize,
    pub polarity: Option<TqPolarity>,
    pub precursor_mz: f64,
    pub product_mz: f64,
    pub time_sec: Vec<f32>,
    pub intensity: Vec<f32>,
}

/// Reader for MRM functions in a TQ MassLynx bundle.
#[derive(Debug, Clone)]
pub struct TqMrmReader {
    pub dir: PathBuf,
    pub functions: Vec<TqMrmFunctionEntry>,
}

impl TqMrmReader {
    /// Discover and validate every MRM function in a bundle.
    ///
    /// A function declared as MRM must have its paired IDX/DAT files. This
    /// reader is the explicit MRM path, so incomplete targeted data is treated
    /// as an error rather than silently dropped.
    pub fn open<P: AsRef<Path>>(dir: P) -> crate::Result<Self> {
        let dir = dir.as_ref().to_path_buf();
        let function_bytes = fs::read(dir.join("_FUNCTNS.INF"))?;
        if function_bytes.len() % FUNCTION_RECORD_SIZE != 0 {
            return Err(crate::Error::Parse(format!(
                "TQ MRM reader: _FUNCTNS.INF size {} is not a multiple of {FUNCTION_RECORD_SIZE}",
                function_bytes.len()
            )));
        }

        let mut functions = Vec::new();
        for (zero_based, record) in function_bytes.chunks_exact(FUNCTION_RECORD_SIZE).enumerate() {
            let index = (zero_based + 1) as u32;
            let descriptor = TqFunctionDescriptor::from_record(record)?;
            if descriptor.kind != TqFunctionKind::Mrm {
                continue;
            }

            let idx_path = dir.join(format!("_FUNC{index:03}.IDX"));
            let dat_path = dir.join(format!("_FUNC{index:03}.DAT"));
            if !idx_path.exists() || !dat_path.exists() {
                return Err(crate::Error::Parse(format!(
                    "TQ MRM reader: function {index} is missing its IDX or DAT file"
                )));
            }

            let idx_bytes = fs::read(&idx_path)?;
            let scan_index = parse_idx22(&idx_bytes)?;
            let dat_size = fs::metadata(&dat_path)?.len();
            let bytes_per_pair = infer_bytes_per_pair(dat_size, &scan_index)?;
            if bytes_per_pair != 4 {
                return Err(crate::Error::Parse(format!(
                    "TQ MRM reader: function {index} uses {bytes_per_pair}-byte records; expected 4-byte packed intensities"
                )));
            }

            validate_mrm_layout(index, &scan_index, bytes_per_pair, dat_size)?;
            let transition_count = consistent_transition_count(index, &scan_index)?;
            let transitions = transitions_from_descriptor(&descriptor, transition_count)?;

            functions.push(TqMrmFunctionEntry {
                index,
                descriptor,
                scan_index,
                dat_path,
                dat_size,
                bytes_per_pair,
                transitions,
            });
        }

        Ok(Self { dir, functions })
    }

    pub fn transition_count(&self) -> usize {
        self.functions
            .iter()
            .map(|function| function.transitions.len())
            .sum()
    }

    /// Decode all targeted functions into one chromatogram per Q1 -> Q3
    /// transition. Retention time is taken from the paired IDX record.
    pub fn chromatograms(&self) -> crate::Result<Vec<TqMrmChromatogram>> {
        let mut output = Vec::with_capacity(self.transition_count());

        for function in &self.functions {
            let dat = fs::read(&function.dat_path)?;
            let mut traces: Vec<Vec<f32>> = function
                .transitions
                .iter()
                .map(|_| Vec::with_capacity(function.cycle_count()))
                .collect();
            let mut times = Vec::with_capacity(function.cycle_count());

            for cycle_index in 0..function.cycle_count() {
                let index_record = &function.scan_index[cycle_index];
                if index_record.pair_count == 0 {
                    continue;
                }

                let range = scan_byte_range(
                    &function.scan_index,
                    cycle_index,
                    function.bytes_per_pair,
                    function.dat_size,
                )?;
                let values = decode_mrm4_cycle(&dat[range])?;
                if values.len() != function.transitions.len() {
                    return Err(crate::Error::Parse(format!(
                        "TQ MRM reader: function {} cycle {} has {} intensity values but {} transitions",
                        function.index,
                        cycle_index,
                        values.len(),
                        function.transitions.len()
                    )));
                }

                times.push(index_record.retention_time_min * 60.0);
                for (trace, value) in traces.iter_mut().zip(values) {
                    trace.push(value);
                }
            }

            for (transition, intensity) in function.transitions.iter().zip(traces) {
                output.push(TqMrmChromatogram {
                    function_index: function.index,
                    transition_index: transition.slot,
                    polarity: function.descriptor.polarity,
                    precursor_mz: transition.precursor_mz,
                    product_mz: transition.product_mz,
                    time_sec: times.clone(),
                    intensity,
                });
            }
        }

        Ok(output)
    }
}

fn validate_mrm_layout(
    function_index: u32,
    scan_index: &[TqIndexRecord],
    bytes_per_pair: usize,
    dat_size: u64,
) -> crate::Result<()> {
    let mut previous_rt: Option<f32> = None;
    for (cycle_index, record) in scan_index.iter().enumerate() {
        if !record.retention_time_min.is_finite() || record.retention_time_min < 0.0 {
            return Err(crate::Error::Parse(format!(
                "TQ MRM reader: function {function_index} cycle {} has invalid retention time {}",
                cycle_index + 1,
                record.retention_time_min
            )));
        }
        if let Some(previous) = previous_rt {
            if record.retention_time_min < previous {
                return Err(crate::Error::Parse(format!(
                    "TQ MRM reader: function {function_index} retention time decreases from {previous} to {} at cycle {}",
                    record.retention_time_min,
                    cycle_index + 1
                )));
            }
        }
        previous_rt = Some(record.retention_time_min);

        let range = scan_byte_range(scan_index, cycle_index, bytes_per_pair, dat_size)?;
        if let Some(next) = scan_index.get(cycle_index + 1) {
            let next_offset = next.dat_offset as usize;
            if next_offset < range.end {
                return Err(crate::Error::Parse(format!(
                    "TQ MRM reader: function {function_index} cycle {} overlaps the next DAT range",
                    cycle_index + 1
                )));
            }
        }
    }
    Ok(())
}

fn consistent_transition_count(
    function_index: u32,
    scan_index: &[TqIndexRecord],
) -> crate::Result<usize> {
    let mut count: Option<u32> = None;
    for record in scan_index.iter().filter(|record| record.pair_count != 0) {
        match count {
            None => count = Some(record.pair_count),
            Some(previous) if previous == record.pair_count => {}
            Some(previous) => {
                return Err(crate::Error::Parse(format!(
                    "TQ MRM reader: function {function_index} changes transition count from {previous} to {} across cycles; scheduled/variable-channel MRM is not yet supported",
                    record.pair_count
                )));
            }
        }
    }

    let count = count.ok_or_else(|| {
        crate::Error::Parse(format!(
            "TQ MRM reader: function {function_index} has no non-empty cycles"
        ))
    })? as usize;

    if count > TRANSITION_SLOTS {
        return Err(crate::Error::Parse(format!(
            "TQ MRM reader: function {function_index} declares {count} transitions but only {TRANSITION_SLOTS} descriptor slots exist"
        )));
    }
    Ok(count)
}

fn transitions_from_descriptor(
    descriptor: &TqFunctionDescriptor,
    count: usize,
) -> crate::Result<Vec<TqMrmTransition>> {
    let mut transitions = Vec::with_capacity(count);
    for slot in 0..count {
        let precursor_mz = descriptor.q1_masses[slot] as f64;
        let product_mz = descriptor.q3_masses[slot] as f64;
        if !precursor_mz.is_finite()
            || !product_mz.is_finite()
            || precursor_mz <= 0.0
            || product_mz <= 0.0
        {
            return Err(crate::Error::Parse(format!(
                "TQ MRM descriptor: transition slot {slot} has invalid Q1/Q3 masses"
            )));
        }
        transitions.push(TqMrmTransition {
            slot,
            precursor_mz,
            product_mz,
        });
    }
    Ok(transitions)
}

/// Decode one acquisition cycle of 4-byte packed MRM intensities.
pub fn decode_mrm4_cycle(cycle_bytes: &[u8]) -> crate::Result<Vec<f32>> {
    if cycle_bytes.len() % 4 != 0 {
        return Err(crate::Error::Parse(format!(
            "TQ MRM 4-byte decoder: cycle size {} is not a multiple of 4",
            cycle_bytes.len()
        )));
    }

    let mut output = Vec::with_capacity(cycle_bytes.len() / 4);
    for record in cycle_bytes.chunks_exact(4) {
        let raw = u32::from_le_bytes([record[0], record[1], record[2], record[3]]);
        let power = (raw >> 22) as i32;
        let base = raw & 0x001f_ffff;
        let intensity = (base as f64 / 1024.0) * 2.0_f64.powi(power - 10);
        output.push(intensity as f32);
    }
    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn packed_mrm_value(power: u32, base: u32) -> [u8; 4] {
        ((power << 22) | (base & 0x001f_ffff)).to_le_bytes()
    }

    fn idx_record(offset: u32, count: u32, rt_min: f32) -> [u8; 22] {
        let mut record = [0_u8; 22];
        record[0..4].copy_from_slice(&offset.to_le_bytes());
        record[4..8].copy_from_slice(&(0x0800_0000_u32 | count).to_le_bytes());
        record[12..16].copy_from_slice(&rt_min.to_le_bytes());
        record
    }

    #[test]
    fn decodes_synthetic_4_byte_intensities() {
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&packed_mrm_value(10, 1024));
        bytes.extend_from_slice(&packed_mrm_value(12, 1024));
        bytes.extend_from_slice(&packed_mrm_value(11, 1536));

        let values = decode_mrm4_cycle(&bytes).unwrap();
        assert_eq!(values, vec![1.0, 4.0, 3.0]);
    }

    #[test]
    fn builds_transition_list_from_synthetic_descriptor() {
        let mut record = [0_u8; FUNCTION_RECORD_SIZE];
        record[0] = 0x09;
        record[0x0a0..0x0a4].copy_from_slice(&300.0_f32.to_le_bytes());
        record[0x0a4..0x0a8].copy_from_slice(&300.0_f32.to_le_bytes());
        record[0x120..0x124].copy_from_slice(&100.0_f32.to_le_bytes());
        record[0x124..0x128].copy_from_slice(&150.0_f32.to_le_bytes());
        let descriptor = TqFunctionDescriptor::from_record(&record).unwrap();

        let transitions = transitions_from_descriptor(&descriptor, 2).unwrap();
        assert_eq!(transitions.len(), 2);
        assert_eq!(transitions[0].precursor_mz, 300.0);
        assert_eq!(transitions[0].product_mz, 100.0);
        assert_eq!(transitions[1].product_mz, 150.0);
    }

    #[test]
    fn ignores_zero_pair_cycles_when_checking_transition_count() {
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&idx_record(0, 2, 0.0));
        bytes.extend_from_slice(&idx_record(8, 0, 0.1));
        bytes.extend_from_slice(&idx_record(8, 2, 0.2));
        let index = parse_idx22(&bytes).unwrap();
        assert_eq!(consistent_transition_count(1, &index).unwrap(), 2);
    }

    #[test]
    fn rejects_decreasing_mrm_retention_time() {
        let index = vec![
            TqIndexRecord {
                dat_offset: 0,
                packed: 2,
                pair_count: 2,
                retention_time_min: 0.2,
            },
            TqIndexRecord {
                dat_offset: 8,
                packed: 2,
                pair_count: 2,
                retention_time_min: 0.1,
            },
        ];
        assert!(validate_mrm_layout(1, &index, 4, 16).is_err());
    }

    #[test]
    fn detects_inconsistent_transition_count() {
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&idx_record(0, 2, 0.0));
        bytes.extend_from_slice(&idx_record(8, 3, 0.1));
        let index = parse_idx22(&bytes).unwrap();
        assert!(consistent_transition_count(1, &index).is_err());
    }
}
