//! Mixed Waters TQ mzML export: broad Q3 spectra plus canonical MRM chromatograms.
//!
//! Broad Q3 scans are emitted through [`crate::tq_mzml::TqQ3Source`], while
//! MRM functions are represented as PSI-MS selected-reaction-monitoring
//! chromatograms (`MS:1001473`) with explicit precursor and product m/z.
//!
//! An optional MRM pseudo-MS2 projection can additionally place sparse
//! Q1-grouped transition spectra into `spectrumList` for compatibility with
//! spectrum-oriented downstream tools. Canonical SRM chromatograms remain
//! present and authoritative in that mode.

use std::io::Write;
use std::path::Path;

use openmassspec_core as msc;
use openmassspec_core::SpectrumSource;

use crate::raw::tq_mrm::{TqMrmChromatogram, TqMrmReader};
use crate::tq_mrm_spectra::{pseudo_ms2_records, TqMrmSpectrumMode};
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
    mrm_spectrum_mode: TqMrmSpectrumMode,
    /// Precomputed only when the optional pseudo-MS2 compatibility projection
    /// is requested. Keeping it here makes conversion errors fail `open`
    /// instead of being silently discarded during iteration, and gives the
    /// mzML writer an exact `spectrumList` count.
    mrm_spectra: Vec<msc::SpectrumRecord>,
}

impl TqMixedSource {
    /// Open both broad-Q3 and canonical MRM views of the same MassLynx bundle.
    ///
    /// This default does not add MRM-derived spectra; MRM remains solely in
    /// canonical SRM chromatograms.
    pub fn open<P: AsRef<Path>>(dir: P, q3_mode: TqQ3MzmlMode) -> crate::Result<Self> {
        Self::open_with_mrm_spectra(dir, q3_mode, TqMrmSpectrumMode::None)
    }

    /// Open the mixed source with an explicit optional MRM spectrum projection.
    ///
    /// [`TqMrmSpectrumMode::PseudoMs2`] adds sparse pseudo-MS2 spectra grouped
    /// by Q1 while still retaining the canonical SRM chromatograms.
    pub fn open_with_mrm_spectra<P: AsRef<Path>>(
        dir: P,
        q3_mode: TqQ3MzmlMode,
        mrm_spectrum_mode: TqMrmSpectrumMode,
    ) -> crate::Result<Self> {
        let dir = dir.as_ref();
        let q3 = TqQ3Source::open(dir, q3_mode)?;
        let mrm = TqMrmReader::open(dir)?;
        let mrm_traces = mrm.chromatograms()?;
        let mut mrm_spectra = match mrm_spectrum_mode {
            TqMrmSpectrumMode::None => Vec::new(),
            // Indices are reassigned after chronological merging with Q3 scans.
            TqMrmSpectrumMode::PseudoMs2 => pseudo_ms2_records(&mrm, 0)?,
        };
        mrm_spectra.sort_by(|left, right| {
            left.retention_time_sec
                .total_cmp(&right.retention_time_sec)
                .then_with(|| left.native_id.cmp(&right.native_id))
        });

        Ok(Self {
            q3,
            mrm_traces,
            mrm_spectrum_mode,
            mrm_spectra,
        })
    }

    pub fn q3_mode(&self) -> TqQ3MzmlMode {
        self.q3.mode()
    }

    pub fn mrm_spectrum_mode(&self) -> TqMrmSpectrumMode {
        self.mrm_spectrum_mode
    }

    pub fn mrm_chromatogram_count(&self) -> usize {
        self.mrm_traces.len()
    }

    pub fn mrm_spectrum_count(&self) -> usize {
        self.mrm_spectra.len()
    }
}

impl SpectrumSource for TqMixedSource {
    fn run_metadata(&self) -> msc::RunMetadata {
        self.q3.run_metadata()
    }

    fn iter_spectra<'s>(&'s mut self) -> Box<dyn Iterator<Item = msc::SpectrumRecord> + 's> {
        if self.mrm_spectra.is_empty() {
            return self.q3.iter_spectra();
        }

        // Clone the optional compatibility stream before borrowing `q3` for
        // the lifetime of its iterator. This keeps the two field borrows
        // unambiguous for the Rust borrow checker and makes the merge logic
        // independent of field evaluation order.
        let mrm_records = self.mrm_spectra.clone();
        let mut q3 = self.q3.iter_spectra().peekable();
        let mut mrm = mrm_records.into_iter().peekable();
        let mut output_index = 0_usize;

        Box::new(std::iter::from_fn(move || {
            let take_q3 = match (q3.peek(), mrm.peek()) {
                (Some(q3_record), Some(mrm_record)) => {
                    q3_record.retention_time_sec <= mrm_record.retention_time_sec
                }
                (Some(_), None) => true,
                (None, Some(_)) => false,
                (None, None) => return None,
            };

            let mut record = if take_q3 { q3.next()? } else { mrm.next()? };
            record.index = output_index;
            record.scan_number = (output_index + 1) as u32;
            output_index += 1;
            Some(record)
        }))
    }

    fn spectrum_count_hint(&self) -> Option<usize> {
        Some(
            self.q3
                .spectrum_count_hint()
                .unwrap_or(0)
                .saturating_add(self.mrm_spectra.len()),
        )
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

/// Write a mixed TQ acquisition to mzML using canonical MRM chromatograms only.
///
/// Q3 semantics are controlled explicitly by `q3_mode`.
pub fn write_tq_mixed_mzml<P: AsRef<Path>, W: Write>(
    dir: P,
    out: &mut W,
    q3_mode: TqQ3MzmlMode,
) -> crate::Result<()> {
    write_tq_mixed_mzml_with_options(dir, out, q3_mode, TqMrmSpectrumMode::None)
}

/// Write a mixed TQ acquisition with explicit Q3 and MRM spectrum projections.
pub fn write_tq_mixed_mzml_with_options<P: AsRef<Path>, W: Write>(
    dir: P,
    out: &mut W,
    q3_mode: TqQ3MzmlMode,
    mrm_spectrum_mode: TqMrmSpectrumMode,
) -> crate::Result<()> {
    let mut source = TqMixedSource::open_with_mrm_spectra(dir, q3_mode, mrm_spectrum_mode)?;
    msc::write_mzml(&mut source, out).map_err(crate::Error::Io)?;
    Ok(())
}

/// Indexed-mzML equivalent of [`write_tq_mixed_mzml`].
pub fn write_tq_mixed_indexed_mzml<P: AsRef<Path>, W: Write>(
    dir: P,
    out: &mut W,
    q3_mode: TqQ3MzmlMode,
) -> crate::Result<()> {
    write_tq_mixed_indexed_mzml_with_options(dir, out, q3_mode, TqMrmSpectrumMode::None)
}

/// Indexed-mzML with explicit Q3 and MRM spectrum projections.
pub fn write_tq_mixed_indexed_mzml_with_options<P: AsRef<Path>, W: Write>(
    dir: P,
    out: &mut W,
    q3_mode: TqQ3MzmlMode,
    mrm_spectrum_mode: TqMrmSpectrumMode,
) -> crate::Result<()> {
    let mut source = TqMixedSource::open_with_mrm_spectra(dir, q3_mode, mrm_spectrum_mode)?;
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
