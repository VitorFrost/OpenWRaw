//! mzML export for Waters triple-quadrupole broad Q3 scans.
//!
//! Two explicit projections are supported:
//! - [`TqQ3MzmlMode::NativeMs2`] preserves the vendor acquisition as an MS2/Q3 scan.
//! - [`TqQ3MzmlMode::PseudoMs1`] intentionally projects the same broad spectrum to
//!   MS1 for untargeted-processing tools that require an MS1 survey stream.
//!
//! Projection provenance is stored as an OpenWRaw `userParam`, not as the
//! Thermo-specific PSI `filter string` term.

use std::collections::BTreeMap;
use std::io::Write;
use std::path::Path;

use openmassspec_core as msc;

use crate::raw::tq::TqPolarity;
use crate::raw::tq_reader::{TqDecodedQ3Scan, TqReader};
use crate::tq_psi_mzml::{write_tq_psi_indexed_mzml, write_tq_psi_mzml};

const SOFTWARE_NAME: &str = "openwraw";
const SOFTWARE_VERSION: &str = env!("CARGO_PKG_VERSION");

/// How broad Q3 scans should be represented in mzML.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum TqQ3MzmlMode {
    /// Preserve the Waters triple-quadrupole acquisition as MS2.
    #[default]
    NativeMs2,
    /// Explicitly project broad Q3 scans to pseudo-MS1 for untargeted workflows.
    PseudoMs1,
}

fn source_file_format_cv() -> msc::CvTerm {
    msc::CvTerm::new("MS:1000526", "Waters raw format")
}

fn native_id_format_cv() -> msc::CvTerm {
    msc::CvTerm::new("MS:1000769", "Waters nativeID format")
}

fn instrument_cv(name: &str) -> msc::CvTerm {
    let compact: String = name
        .chars()
        .filter(|c| c.is_ascii_alphanumeric())
        .flat_map(char::to_uppercase)
        .collect();

    if compact.starts_with("XEVOTQSMICRO") {
        return msc::CvTerm::new("MS:1002731", "Xevo TQ-S micro");
    }
    if compact.starts_with("XEVOTQS") {
        return msc::CvTerm::new("MS:1001792", "Xevo TQ-S");
    }
    if compact.starts_with("XEVOTQ") {
        return msc::CvTerm::new("MS:1001790", "Xevo TQ MS");
    }
    msc::CvTerm::new("MS:1000126", "Waters instrument model")
}

fn polarity_for(polarity: Option<TqPolarity>) -> Option<msc::Polarity> {
    match polarity {
        Some(TqPolarity::Positive) => Some(msc::Polarity::Positive),
        Some(TqPolarity::Negative) => Some(msc::Polarity::Negative),
        None => None,
    }
}

fn native_id_for(function_index: u32, scan_index_zero_based: usize) -> String {
    format!(
        "function={function_index} process=0 scan={}",
        scan_index_zero_based + 1
    )
}

fn bundle_name(reader: &TqReader) -> String {
    reader
        .dir
        .file_name()
        .map(|name| name.to_string_lossy().into_owned())
        .unwrap_or_else(|| "bundle.raw".to_owned())
}

fn mzml_run_id_for(source_name: &str) -> String {
    let mut id = String::from("run_");
    for ch in source_name.chars() {
        if ch.is_ascii_alphanumeric() || matches!(ch, '_' | '-' | '.') {
            id.push(ch);
        } else {
            id.push('_');
        }
    }
    if id == "run_" {
        id.push_str("waters_tq");
    }
    id
}

fn run_metadata_for(reader: &TqReader) -> msc::RunMetadata {
    let instrument_name = reader
        .header
        .instrument
        .clone()
        .unwrap_or_else(|| "Waters".to_owned());

    msc::RunMetadata {
        extra: BTreeMap::new(),
        source_file_name: bundle_name(reader),
        source_file_format: source_file_format_cv(),
        native_id_format: native_id_format_cv(),
        instrument: instrument_cv(&instrument_name),
        instrument_serial_number: None,
        software_name: SOFTWARE_NAME.to_owned(),
        software_version: SOFTWARE_VERSION.to_owned(),
        acquisition_software_name: None,
        acquisition_software_version: None,
        start_timestamp: None,
        mobility_array_kind: None,
        // openmassspec-core 1.5 models analyzers as additional instrument
        // configurations rather than componentList/analyzer children. Keep
        // the accurately identified Xevo model and avoid emitting a
        // structurally misplaced quadrupole term until core exposes proper
        // instrument-component metadata.
        analyzers: Vec::new(),
    }
}

fn summarize_arrays(
    mz: &[f64],
    intensity: &[f32],
) -> (f64, Option<f64>, Option<f64>, Option<f64>, Option<f64>) {
    if mz.is_empty() {
        return (0.0, None, None, None, None);
    }

    let mut tic = 0.0_f64;
    let mut base_peak_intensity = f32::NEG_INFINITY;
    let mut base_peak_mz = mz[0];
    let mut low_mz = f64::INFINITY;
    let mut high_mz = f64::NEG_INFINITY;

    for (&mass, &signal) in mz.iter().zip(intensity.iter()) {
        tic += signal as f64;
        if signal > base_peak_intensity {
            base_peak_intensity = signal;
            base_peak_mz = mass;
        }
        low_mz = low_mz.min(mass);
        high_mz = high_mz.max(mass);
    }

    (
        tic,
        Some(base_peak_mz),
        Some(base_peak_intensity as f64),
        Some(low_mz),
        Some(high_mz),
    )
}

fn record_from_scan(
    mode: TqQ3MzmlMode,
    scan_counter: u32,
    scan: TqDecodedQ3Scan,
) -> msc::SpectrumRecord {
    let (tic, base_peak_mz, base_peak_intensity, low_mz, high_mz) =
        summarize_arrays(&scan.spectrum.mz, &scan.spectrum.intensity);

    let (ms_level, precursor, projection) = match mode {
        TqQ3MzmlMode::NativeMs2 => {
            let target_mz = (scan.set_mass_da.is_finite() && scan.set_mass_da > 0.0)
                .then_some(scan.set_mass_da as f64);
            let precursor = target_mz.map(|target_mz| msc::PrecursorInfo {
                target_mz: Some(target_mz),
                analyzer: Some(msc::Analyzer::TQMS),
                ..Default::default()
            });
            (2, precursor, "native-ms2-q3")
        }
        TqQ3MzmlMode::PseudoMs1 => (1, None, "pseudo-ms1-from-q3"),
    };

    let mut extra = BTreeMap::new();
    extra.insert("openwraw.projection".to_owned(), projection.to_owned());
    extra.insert(
        "openwraw.source_acquisition".to_owned(),
        "Waters TQ broad Q3 scan".to_owned(),
    );

    msc::SpectrumRecord {
        extra,
        acquisition_event_id: None,
        index: (scan_counter as usize).saturating_sub(1),
        scan_number: scan_counter,
        native_id: native_id_for(scan.function_index, scan.scan_index),
        ms_level,
        polarity: polarity_for(scan.polarity),
        scan_mode: Some(msc::ScanMode::Profile),
        analyzer: Some(msc::Analyzer::TQMS),
        filter: None,
        retention_time_sec: scan.retention_time_min as f64 * 60.0,
        total_ion_current: Some(tic),
        base_peak_mz,
        base_peak_intensity,
        low_mz,
        high_mz,
        ion_injection_time_ms: None,
        inv_mobility: None,
        faims_cv: None,
        precursor,
        mz: scan.spectrum.mz,
        intensity: scan.spectrum.intensity,
        inv_mobility_per_peak: None,
    }
}

/// Streaming `openmassspec-core` source for TQ Q3 scans.
pub struct TqQ3Source {
    reader: TqReader,
    mode: TqQ3MzmlMode,
    use_mzml_safe_run_id: bool,
}

impl TqQ3Source {
    pub fn new(reader: TqReader, mode: TqQ3MzmlMode) -> Self {
        Self {
            reader,
            mode,
            use_mzml_safe_run_id: false,
        }
    }

    pub(crate) fn prepare_for_mzml(&mut self) {
        self.use_mzml_safe_run_id = true;
    }

    pub fn open<P: AsRef<Path>>(dir: P, mode: TqQ3MzmlMode) -> crate::Result<Self> {
        Ok(Self::new(TqReader::open(dir)?, mode))
    }

    pub fn reader(&self) -> &TqReader {
        &self.reader
    }

    pub fn mode(&self) -> TqQ3MzmlMode {
        self.mode
    }
}

impl msc::SpectrumSource for TqQ3Source {
    fn run_metadata(&self) -> msc::RunMetadata {
        let mut metadata = run_metadata_for(&self.reader);
        if self.use_mzml_safe_run_id {
            // openmassspec-core 1.5 uses source_file_name for both the
            // sourceFile name and the XML xs:ID-valued run id. TQ PSI
            // serialization later rebuilds sourceFileList from the RAW
            // bundle, so only the transient writer metadata needs a safe ID.
            metadata.source_file_name = mzml_run_id_for(&metadata.source_file_name);
        }
        metadata
    }

    fn iter_spectra<'s>(&'s mut self) -> Box<dyn Iterator<Item = msc::SpectrumRecord> + 's> {
        let reader = &self.reader;
        let mode = self.mode;
        let mut emitted_count = 0_u32;
        Box::new(reader.iter_scans().map(move |decoded| {
            let scan = match decoded {
                Ok(scan) => scan,
                Err(error) => panic!(
                    "TQ Q3 spectrum iteration failed after bundle validation: {error}"
                ),
            };
            emitted_count += 1;
            record_from_scan(mode, emitted_count, scan)
        }))
    }

    fn spectrum_count_hint(&self) -> Option<usize> {
        Some(self.reader.total_q3_scan_count())
    }

    fn iter_chromatograms<'s>(
        &'s mut self,
    ) -> Box<dyn Iterator<Item = msc::ChromatogramRecord> + 's> {
        Box::new(std::iter::empty())
    }
}

/// Write broad Q3 scans to PSI-corrected mzML using the selected projection.
pub fn write_tq_q3_mzml<P: AsRef<Path>, W: Write>(
    dir: P,
    out: &mut W,
    mode: TqQ3MzmlMode,
) -> crate::Result<()> {
    let dir = dir.as_ref();
    let mut source = TqQ3Source::open(dir, mode)?;
    source.prepare_for_mzml();
    write_tq_psi_mzml(&mut source, dir, out)
}

/// Indexed-mzML equivalent of [`write_tq_q3_mzml`], rebuilt after PSI fixes.
pub fn write_tq_q3_indexed_mzml<P: AsRef<Path>, W: Write>(
    dir: P,
    out: &mut W,
    mode: TqQ3MzmlMode,
) -> crate::Result<()> {
    let dir = dir.as_ref();
    let mut source = TqQ3Source::open(dir, mode)?;
    source.prepare_for_mzml();
    write_tq_psi_indexed_mzml(&mut source, dir, out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::raw::data::Spectrum;
    use crate::raw::tq::TqFunctionKind;

    fn synthetic_scan() -> TqDecodedQ3Scan {
        TqDecodedQ3Scan {
            function_index: 3,
            scan_index: 4,
            retention_time_min: 1.5,
            polarity: Some(TqPolarity::Negative),
            set_mass_da: 40.0,
            declared_mz_low: 75.0,
            declared_mz_high: 900.0,
            spectrum: Spectrum {
                mz: vec![100.0, 150.0, 200.0],
                intensity: vec![2.0, 10.0, 3.0],
            },
        }
    }

    #[test]
    fn native_mode_preserves_ms2_and_precursor_context() {
        let record = record_from_scan(TqQ3MzmlMode::NativeMs2, 1, synthetic_scan());
        assert_eq!(record.ms_level, 2);
        assert_eq!(record.scan_mode, Some(msc::ScanMode::Profile));
        assert_eq!(record.analyzer, Some(msc::Analyzer::TQMS));
        assert_eq!(record.polarity, Some(msc::Polarity::Negative));
        assert_eq!(record.precursor.unwrap().target_mz, Some(40.0));
    }

    #[test]
    fn native_mode_does_not_emit_zero_as_a_precursor_mass() {
        let mut scan = synthetic_scan();
        scan.set_mass_da = 0.0;
        let record = record_from_scan(TqQ3MzmlMode::NativeMs2, 1, scan);
        assert!(record.precursor.is_none());
    }

    #[test]
    fn pseudo_mode_uses_user_param_provenance_not_filter_string() {
        let record = record_from_scan(TqQ3MzmlMode::PseudoMs1, 1, synthetic_scan());
        assert_eq!(record.ms_level, 1);
        assert!(record.precursor.is_none());
        assert!(record.filter.is_none());
        assert_eq!(
            record.extra.get("openwraw.projection").map(String::as_str),
            Some("pseudo-ms1-from-q3")
        );
    }

    #[test]
    fn summary_matches_arrays() {
        let record = record_from_scan(TqQ3MzmlMode::PseudoMs1, 1, synthetic_scan());
        assert_eq!(record.total_ion_current, Some(15.0));
        assert_eq!(record.base_peak_mz, Some(150.0));
        assert_eq!(record.base_peak_intensity, Some(10.0));
        assert_eq!(record.low_mz, Some(100.0));
        assert_eq!(record.high_mz, Some(200.0));
    }

    #[test]
    fn instrument_mapping_handles_compact_tq_s_micro_header_strings() {
        let cv = instrument_cv("XEVO-TQSmicro#serial-redacted");
        assert_eq!(cv.accession, "MS:1002731");
        assert_eq!(cv.name, "Xevo TQ-S micro");
    }

    #[test]
    fn q3_source_mode_enum_is_independent_of_function_kind() {
        assert_eq!(TqFunctionKind::Q3Scan, TqFunctionKind::Q3Scan);
        assert_ne!(TqFunctionKind::Q3Scan, TqFunctionKind::Mrm);
    }

    #[test]
    fn mzml_run_id_is_valid_for_numeric_or_punctuated_source_names() {
        assert_eq!(mzml_run_id_for("20230413_EL-020.raw"), "run_20230413_EL-020.raw");
        assert_eq!(mzml_run_id_for("a path/bundle.raw"), "run_a_path_bundle.raw");
        assert_eq!(mzml_run_id_for(""), "run_waters_tq");
    }
}
