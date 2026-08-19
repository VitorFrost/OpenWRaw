//! Mixed Waters TQ mzML export: broad Q3 spectra plus canonical MRM chromatograms.
//!
//! Broad Q3 scans are emitted through [`crate::tq_mzml::TqQ3Source`], while
//! MRM functions are represented as PSI-MS selected-reaction-monitoring
//! chromatograms (`MS:1001473`) with explicit precursor and product m/z.

use std::io::Write;
use std::path::Path;

use openmassspec_core as msc;
use openmassspec_core::SpectrumSource;

use crate::raw::tq_mrm::{TqMrmChromatogram, TqMrmReader};
use crate::tq_mzml::{TqQ3MzmlMode, TqQ3Source};

fn mrm_record(index: usize, trace: TqMrmChromatogram) -> msc::ChromatogramRecord {
    msc::ChromatogramRecord {
        index,
        id: format!(
            "function={} transition={}",
            trace.function_index,
            trace.transition_index + 1
        ),
        chromatogram_type: Some(msc::CvTerm::new(
            "MS:1001473",
            "selected reaction monitoring chromatogram",
        )),
        precursor_mz: Some(trace.precursor_mz),
        product_mz: Some(trace.product_mz),
        time_sec: trace.time_sec,
        intensity: trace.intensity,
    }
}

/// Spectrum/chromatogram source for a mixed TQ acquisition.
pub struct TqMixedSource {
    q3: TqQ3Source,
    mrm_traces: Vec<TqMrmChromatogram>,
}

impl TqMixedSource {
    /// Open both broad-Q3 and MRM views of the same MassLynx bundle.
    ///
    /// MRM decoding errors are surfaced here rather than silently dropping
    /// targeted channels from the resulting mzML.
    pub fn open<P: AsRef<Path>>(dir: P, q3_mode: TqQ3MzmlMode) -> crate::Result<Self> {
        let dir = dir.as_ref();
        let q3 = TqQ3Source::open(dir, q3_mode)?;
        let mrm = TqMrmReader::open(dir)?;
        let mrm_traces = mrm.chromatograms()?;
        Ok(Self { q3, mrm_traces })
    }

    pub fn q3_mode(&self) -> TqQ3MzmlMode {
        self.q3.mode()
    }

    pub fn mrm_chromatogram_count(&self) -> usize {
        self.mrm_traces.len()
    }
}

impl SpectrumSource for TqMixedSource {
    fn run_metadata(&self) -> msc::RunMetadata {
        self.q3.run_metadata()
    }

    fn iter_spectra<'s>(&'s mut self) -> Box<dyn Iterator<Item = msc::SpectrumRecord> + 's> {
        self.q3.iter_spectra()
    }

    fn spectrum_count_hint(&self) -> Option<usize> {
        self.q3.spectrum_count_hint()
    }

    fn iter_chromatograms<'s>(
        &'s mut self,
    ) -> Box<dyn Iterator<Item = msc::ChromatogramRecord> + 's> {
        let records: Vec<msc::ChromatogramRecord> = self
            .mrm_traces
            .clone()
            .into_iter()
            .enumerate()
            .map(|(index, trace)| mrm_record(index, trace))
            .collect();
        Box::new(records.into_iter())
    }
}

/// Write a mixed TQ acquisition to mzML.
///
/// Q3 semantics are controlled explicitly by `q3_mode`; MRM data always
/// remains canonical SRM chromatograms.
pub fn write_tq_mixed_mzml<P: AsRef<Path>, W: Write>(
    dir: P,
    out: &mut W,
    q3_mode: TqQ3MzmlMode,
) -> crate::Result<()> {
    let mut source = TqMixedSource::open(dir, q3_mode)?;
    msc::write_mzml(&mut source, out).map_err(crate::Error::Io)?;
    Ok(())
}

/// Indexed-mzML equivalent of [`write_tq_mixed_mzml`].
pub fn write_tq_mixed_indexed_mzml<P: AsRef<Path>, W: Write>(
    dir: P,
    out: &mut W,
    q3_mode: TqQ3MzmlMode,
) -> crate::Result<()> {
    let mut source = TqMixedSource::open(dir, q3_mode)?;
    msc::write_indexed_mzml(&mut source, out).map_err(crate::Error::Io)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::raw::tq::TqPolarity;

    #[test]
    fn mrm_record_uses_canonical_srm_chromatogram_term() {
        let record = mrm_record(
            2,
            TqMrmChromatogram {
                function_index: 4,
                transition_index: 1,
                polarity: Some(TqPolarity::Positive),
                precursor_mz: 300.0,
                product_mz: 150.0,
                time_sec: vec![0.0, 1.0],
                intensity: vec![10.0, 20.0],
            },
        );

        assert_eq!(record.index, 2);
        assert_eq!(record.precursor_mz, Some(300.0));
        assert_eq!(record.product_mz, Some(150.0));
        let cv = record.chromatogram_type.unwrap();
        assert_eq!(cv.accession, "MS:1001473");
        assert_eq!(cv.name, "selected reaction monitoring chromatogram");
    }
}
