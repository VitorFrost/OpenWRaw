//! High-level reader for broad Q3 scans in Waters triple-quadrupole RAW bundles.
//!
//! This reader is intentionally separate from the existing TOF [`crate::reader::Reader`]
//! while TQ support is being validated. It does not require `_extern.inf` TOF geometry
//! and it does not silently reinterpret MRM functions as ordinary spectra.

use std::fs;
use std::io::{Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};

use crate::raw::data::Spectrum;
use crate::raw::header::{FunctionCal, Header};
use crate::raw::tq::{
    decode_direct6, infer_bytes_per_pair, parse_idx22, scan_byte_range, TqFunctionDescriptor,
    TqFunctionKind, TqIndexRecord, TqPolarity, FUNCTION_RECORD_SIZE,
};

/// One broad Q3-scan function ready for random-access decoding.
#[derive(Debug, Clone)]
pub struct TqQ3FunctionEntry {
    pub index: u32,
    pub descriptor: TqFunctionDescriptor,
    pub scan_index: Vec<TqIndexRecord>,
    pub dat_path: PathBuf,
    pub dat_size: u64,
    pub bytes_per_pair: usize,
    pub calibration: FunctionCal,
}

impl TqQ3FunctionEntry {
    pub fn scan_count(&self) -> usize {
        self.scan_index.len()
    }
}

/// One decoded broad Q3 scan.
#[derive(Debug, Clone)]
pub struct TqDecodedQ3Scan {
    pub function_index: u32,
    pub scan_index: usize,
    pub retention_time_min: f32,
    pub polarity: Option<TqPolarity>,
    /// Function-level set-mass / precursor-context metadata.
    pub set_mass_da: f32,
    pub declared_mz_low: f32,
    pub declared_mz_high: f32,
    pub spectrum: Spectrum,
}

/// Reader for the broad Q3-scan subset of a TQ MassLynx bundle.
///
/// MRM functions are counted and exposed through [`Self::mrm_function_count`],
/// but are not converted to spectra by this reader. Canonical MRM support needs
/// chromatogram semantics and is implemented separately from broad-scan decoding.
#[derive(Debug, Clone)]
pub struct TqReader {
    pub dir: PathBuf,
    pub header: Header,
    pub q3_functions: Vec<TqQ3FunctionEntry>,
    pub mrm_function_count: usize,
    pub unknown_function_count: usize,
}

impl TqReader {
    /// Open a Waters `.raw` directory and discover TQ Q3-scan functions.
    ///
    /// Unlike the TOF reader, this path intentionally does not require
    /// `_extern.inf` flight-path geometry because direct low-resolution packed
    /// data stores m/z rather than TOF-bin coordinates.
    pub fn open<P: AsRef<Path>>(dir: P) -> crate::Result<Self> {
        let dir = dir.as_ref().to_path_buf();
        let header = Header::from_path(&dir.join("_HEADER.TXT"))?;
        let function_bytes = fs::read(dir.join("_FUNCTNS.INF"))?;
        if function_bytes.len() % FUNCTION_RECORD_SIZE != 0 {
            return Err(crate::Error::Parse(format!(
                "TQ reader: _FUNCTNS.INF size {} is not a multiple of {FUNCTION_RECORD_SIZE}",
                function_bytes.len()
            )));
        }

        let mut q3_functions = Vec::new();
        let mut mrm_function_count = 0_usize;
        let mut unknown_function_count = 0_usize;

        for (zero_based, record) in function_bytes.chunks_exact(FUNCTION_RECORD_SIZE).enumerate() {
            let index = (zero_based + 1) as u32;
            let descriptor = TqFunctionDescriptor::from_record(record)?;
            match descriptor.kind {
                TqFunctionKind::Mrm => {
                    mrm_function_count += 1;
                }
                TqFunctionKind::Q3Scan => {
                    let idx_path = dir.join(format!("_FUNC{index:03}.IDX"));
                    let dat_path = dir.join(format!("_FUNC{index:03}.DAT"));
                    if !idx_path.exists() || !dat_path.exists() {
                        return Err(crate::Error::Parse(format!(
                            "TQ reader: Q3 function {index} is missing its IDX or DAT file"
                        )));
                    }

                    let idx_bytes = fs::read(&idx_path)?;
                    let scan_index = parse_idx22(&idx_bytes)?;
                    let dat_size = fs::metadata(&dat_path)?.len();
                    let bytes_per_pair = infer_bytes_per_pair(dat_size, &scan_index)?;
                    if bytes_per_pair != 6 {
                        return Err(crate::Error::Parse(format!(
                            "TQ reader: Q3 function {index} uses {bytes_per_pair}-byte records; expected direct 6-byte low-resolution encoding"
                        )));
                    }
                    validate_q3_layout(index, &scan_index, bytes_per_pair, dat_size)?;

                    let calibration = header
                        .cal_functions
                        .get(&index)
                        .cloned()
                        .unwrap_or_default();

                    q3_functions.push(TqQ3FunctionEntry {
                        index,
                        descriptor,
                        scan_index,
                        dat_path,
                        dat_size,
                        bytes_per_pair,
                        calibration,
                    });
                }
                TqFunctionKind::Unknown(_) => {
                    unknown_function_count += 1;
                }
            }
        }

        Ok(Self {
            dir,
            header,
            q3_functions,
            mrm_function_count,
            unknown_function_count,
        })
    }

    pub fn total_q3_scan_count(&self) -> usize {
        self.q3_functions
            .iter()
            .map(TqQ3FunctionEntry::scan_count)
            .sum()
    }

    /// Decode one scan from one Q3 function.
    pub fn decode_scan(
        &self,
        function_index: u32,
        scan_index: usize,
    ) -> crate::Result<TqDecodedQ3Scan> {
        let function = self
            .q3_functions
            .iter()
            .find(|entry| entry.index == function_index)
            .ok_or_else(|| {
                crate::Error::Parse(format!(
                    "TQ reader: Q3 function {function_index} is not present"
                ))
            })?;

        let range = scan_byte_range(
            &function.scan_index,
            scan_index,
            function.bytes_per_pair,
            function.dat_size,
        )?;
        let bytes = read_range(&function.dat_path, range.start as u64, range.len())?;
        let spectrum = decode_direct6(&bytes, Some(&function.calibration))?;
        let idx = &function.scan_index[scan_index];

        Ok(TqDecodedQ3Scan {
            function_index,
            scan_index,
            retention_time_min: idx.retention_time_min,
            polarity: function.descriptor.polarity,
            set_mass_da: function.descriptor.set_mass_da,
            declared_mz_low: function.descriptor.mz_low,
            declared_mz_high: function.descriptor.mz_high,
            spectrum,
        })
    }

    /// Iterate all Q3 scans in acquisition order across functions.
    ///
    /// MassLynx stores each function in a separate `_FUNCnnn.*` stream even
    /// when functions were interleaved during acquisition. Ordering only by
    /// function would therefore turn a polarity-switching run into separate
    /// full-duration blocks. Retention time reconstructs the acquisition
    /// chronology while native IDs still retain the original function/scan.
    pub fn iter_scans(&self) -> impl Iterator<Item = crate::Result<TqDecodedQ3Scan>> + '_ {
        let mut plan: Vec<(f32, u32, usize)> = self
            .q3_functions
            .iter()
            .flat_map(|function| {
                function
                    .scan_index
                    .iter()
                    .enumerate()
                    .map(move |(scan, record)| {
                        (record.retention_time_min, function.index, scan)
                    })
            })
            .collect();

        plan.sort_by(|left, right| {
            left.0
                .total_cmp(&right.0)
                .then_with(|| left.1.cmp(&right.1))
                .then_with(|| left.2.cmp(&right.2))
        });

        plan.into_iter()
            .map(move |(_, function, scan)| self.decode_scan(function, scan))
    }
}

/// Validate all static IDX -> DAT relationships while opening the bundle.
///
/// `decode_direct6` cannot fail for a correctly-sized six-byte range, so this
/// catches deterministic truncation/overlap before `SpectrumSource` advertises
/// its spectrum count. Runtime I/O failure or a file being modified after open
/// remains a separate environmental error.
fn validate_q3_layout(
    function_index: u32,
    scan_index: &[TqIndexRecord],
    bytes_per_pair: usize,
    dat_size: u64,
) -> crate::Result<()> {
    for (scan_number, record) in scan_index.iter().enumerate() {
        if !record.retention_time_min.is_finite() || record.retention_time_min < 0.0 {
            return Err(crate::Error::Parse(format!(
                "TQ reader: function {function_index} scan {} has invalid retention time {}",
                scan_number + 1,
                record.retention_time_min
            )));
        }

        let range = scan_byte_range(scan_index, scan_number, bytes_per_pair, dat_size)?;
        if let Some(next) = scan_index.get(scan_number + 1) {
            let next_offset = next.dat_offset as usize;
            if next_offset < range.end {
                return Err(crate::Error::Parse(format!(
                    "TQ reader: function {function_index} scan {} overlaps the next DAT range (end {}, next offset {next_offset})",
                    scan_number + 1,
                    range.end
                )));
            }
        }
    }
    Ok(())
}

fn read_range(path: &Path, offset: u64, length: usize) -> crate::Result<Vec<u8>> {
    let mut file = fs::File::open(path)?;
    file.seek(SeekFrom::Start(offset))?;
    let mut bytes = vec![0_u8; length];
    file.read_exact(&mut bytes)?;
    Ok(bytes)
}

#[cfg(test)]
mod tests {
    use super::*;

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

    fn make_idx_record(offset: u32, pair_count: u32, rt_min: f32) -> [u8; 22] {
        let mut record = [0_u8; 22];
        record[0..4].copy_from_slice(&offset.to_le_bytes());
        let packed = 0x1800_0000_u32 | pair_count;
        record[4..8].copy_from_slice(&packed.to_le_bytes());
        record[12..16].copy_from_slice(&rt_min.to_le_bytes());
        record
    }

    fn synthetic_bundle_dir() -> PathBuf {
        std::env::temp_dir().join(format!("openwraw-tq-reader-{}", std::process::id()))
    }

    #[test]
    fn opens_and_decodes_synthetic_q3_bundle_without_extern_inf() {
        let dir = synthetic_bundle_dir();
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();

        let header = "$$ Version: 01.00\r\n\
$$ Instrument: SYNTHETIC-TQ\r\n\
$$ Cal Function 1: 0.0,1.0,T0\r\n";
        fs::write(dir.join("_HEADER.TXT"), header).unwrap();

        let mut function = [0_u8; FUNCTION_RECORD_SIZE];
        function[0] = 0x0b;
        function[0x018..0x01c].copy_from_slice(&40.0_f32.to_le_bytes());
        function[0x020..0x024].copy_from_slice(&0.5_f32.to_le_bytes());
        function[0x0a0..0x0a4].copy_from_slice(&75.0_f32.to_le_bytes());
        function[0x120..0x124].copy_from_slice(&900.0_f32.to_le_bytes());
        fs::write(dir.join("_FUNCTNS.INF"), function).unwrap();

        let idx = make_idx_record(0, 2, 0.25);
        fs::write(dir.join("_FUNC001.IDX"), idx).unwrap();

        let mut dat = Vec::new();
        dat.extend_from_slice(&direct6_record(100, 23, 10, 0));
        dat.extend_from_slice(&direct6_record(201, 22, 5, 1));
        fs::write(dir.join("_FUNC001.DAT"), dat).unwrap();

        let reader = TqReader::open(&dir).unwrap();
        assert_eq!(reader.q3_functions.len(), 1);
        assert_eq!(reader.mrm_function_count, 0);
        assert_eq!(reader.total_q3_scan_count(), 1);

        let scan = reader.decode_scan(1, 0).unwrap();
        assert_eq!(scan.polarity, Some(TqPolarity::Positive));
        assert!((scan.retention_time_min - 0.25).abs() < f32::EPSILON);
        assert_eq!(scan.spectrum.mz.len(), 2);
        assert!((scan.spectrum.mz[0] - 100.0).abs() < 1e-12);
        assert!((scan.spectrum.mz[1] - 100.5).abs() < 1e-12);
        assert_eq!(scan.spectrum.intensity, vec![10.0, 20.0]);

        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn rejects_overlapping_q3_scan_ranges() {
        let index = vec![
            TqIndexRecord {
                dat_offset: 0,
                packed: 2,
                pair_count: 2,
                retention_time_min: 0.0,
            },
            TqIndexRecord {
                dat_offset: 6,
                packed: 2,
                pair_count: 2,
                retention_time_min: 0.1,
            },
        ];
        assert!(validate_q3_layout(1, &index, 6, 18).is_err());
    }

    #[test]
    fn counts_mrm_functions_without_reinterpreting_them_as_q3_scans() {
        let dir = std::env::temp_dir().join(format!(
            "openwraw-tq-reader-mrm-count-{}",
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        fs::write(dir.join("_HEADER.TXT"), "$$ Version: 01.00\r\n").unwrap();

        let mut functions = vec![0_u8; FUNCTION_RECORD_SIZE * 2];
        functions[0] = 0x09;
        functions[FUNCTION_RECORD_SIZE] = 0x29;
        fs::write(dir.join("_FUNCTNS.INF"), functions).unwrap();

        let reader = TqReader::open(&dir).unwrap();
        assert_eq!(reader.mrm_function_count, 2);
        assert!(reader.q3_functions.is_empty());
        assert_eq!(reader.total_q3_scan_count(), 0);

        let _ = fs::remove_dir_all(&dir);
    }
}
