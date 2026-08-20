//! Optional pseudo-MS2 projection for Waters TQ MRM functions.
//!
//! Canonical MRM data remains represented as SRM chromatograms. This module
//! exists only for downstream tools that operate primarily on spectrum lists.
//! Each acquisition cycle is grouped by Q1 precursor mass and converted into
//! a sparse MS2 spectrum whose m/z axis is the set of Q3 product masses.

use std::fs;

use openmassspec_core as msc;

use crate::raw::tq::{scan_byte_range, TqPolarity};
use crate::raw::tq_mrm::{decode_mrm4_cycle, TqMrmFunctionEntry, TqMrmReader, TqMrmTransition};

/// Optional MRM representation in the mzML spectrum list.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum TqMrmSpectrumMode {
    /// Do not add MRM-derived spectra. Canonical SRM chromatograms are still emitted.
    #[default]
    None,
    /// Add sparse pseudo-MS2 spectra grouped by Q1 precursor mass.
    PseudoMs2,
}

#[derive(Debug, Clone)]
struct TransitionGroup {
    precursor_mz: f64,
    transition_indices: Vec<usize>,
}

fn transition_groups(transitions: &[TqMrmTransition]) -> Vec<TransitionGroup> {
    let mut groups: Vec<TransitionGroup> = Vec::new();
    for (transition_index, transition) in transitions.iter().enumerate() {
        if let Some(group) = groups
            .iter_mut()
            .find(|group| (group.precursor_mz - transition.precursor_mz).abs() < 1e-6)
        {
            group.transition_indices.push(transition_index);
        } else {
            groups.push(TransitionGroup {
                precursor_mz: transition.precursor_mz,
                transition_indices: vec![transition_index],
            });
        }
    }
    groups
}

fn polarity_for(polarity: Option<TqPolarity>) -> Option<msc::Polarity> {
    match polarity {
        Some(TqPolarity::Positive) => Some(msc::Polarity::Positive),
        Some(TqPolarity::Negative) => Some(msc::Polarity::Negative),
        None => None,
    }
}

fn summarize(
    mz: &[f64],
    intensity: &[f32],
) -> (f64, Option<f64>, Option<f64>, Option<f64>, Option<f64>) {
    if mz.is_empty() {
        return (0.0, None, None, None, None);
    }

    let mut tic = 0.0_f64;
    let mut base_peak_mz = mz[0];
    let mut base_peak_intensity = f32::NEG_INFINITY;
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

fn sorted_group_arrays(
    function: &TqMrmFunctionEntry,
    values: &[f32],
    group: &TransitionGroup,
) -> (Vec<f64>, Vec<f32>) {
    let mut pairs: Vec<(f64, f32)> = group
        .transition_indices
        .iter()
        .map(|&transition_index| {
            (
                function.transitions[transition_index].product_mz,
                values[transition_index],
            )
        })
        .collect();

    // MRM transition order is method-defined and need not be ascending in Q3.
    // A spectrum representation is easier and safer for downstream consumers
    // when m/z is ordered, so sort pairs together and preserve signal pairing.
    pairs.sort_by(|a, b| {
        a.0.partial_cmp(&b.0)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    pairs.into_iter().unzip()
}

fn function_pseudo_ms2(
    function: &TqMrmFunctionEntry,
    index_offset: usize,
) -> crate::Result<Vec<msc::SpectrumRecord>> {
    let dat = fs::read(&function.dat_path)?;
    let groups = transition_groups(&function.transitions);
    let mut output = Vec::with_capacity(function.cycle_count() * groups.len());

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
                "TQ MRM pseudo-MS2: function {} cycle {} has {} values for {} transitions",
                function.index,
                cycle_index,
                values.len(),
                function.transitions.len()
            )));
        }

        for (group_index, group) in groups.iter().enumerate() {
            let (mz, intensity) = sorted_group_arrays(function, &values, group);
            let (tic, base_peak_mz, base_peak_intensity, low_mz, high_mz) =
                summarize(&mz, &intensity);
            let index = index_offset + output.len();
            output.push(msc::SpectrumRecord {
                extra: ::std::collections::BTreeMap::new(),
                acquisition_event_id: None,
                index,
                scan_number: (index + 1) as u32,
                // Keep the declared Waters nativeID token structure. The
                // process number is used only to distinguish synthetic Q1
                // groups from the same physical acquisition cycle.
                native_id: format!(
                    "function={} process={} scan={}",
                    function.index,
                    group_index + 1,
                    cycle_index + 1
                ),
                ms_level: 2,
                polarity: polarity_for(function.descriptor.polarity),
                scan_mode: Some(msc::ScanMode::Centroid),
                analyzer: Some(msc::Analyzer::TQMS),
                filter: Some(
                    "OpenWRaw pseudo-MS2 projection of Waters TQ MRM transitions".to_owned(),
                ),
                retention_time_sec: index_record.retention_time_min as f64 * 60.0,
                total_ion_current: Some(tic),
                base_peak_mz,
                base_peak_intensity,
                low_mz,
                high_mz,
                ion_injection_time_ms: None,
                inv_mobility: None,
                faims_cv: None,
                precursor: Some(msc::PrecursorInfo {
                    target_mz: Some(group.precursor_mz),
                    selected_mz: Some(group.precursor_mz),
                    analyzer: Some(msc::Analyzer::TQMS),
                    ..Default::default()
                }),
                mz,
                intensity,
                inv_mobility_per_peak: None,
            });
        }
    }

    Ok(output)
}

/// Convert all MRM functions into optional sparse pseudo-MS2 spectra.
///
/// `index_offset` is used when these records are appended after another
/// spectrum stream so global mzML spectrum indices remain contiguous.
pub fn pseudo_ms2_records(
    reader: &TqMrmReader,
    index_offset: usize,
) -> crate::Result<Vec<msc::SpectrumRecord>> {
    let mut output = Vec::new();
    for function in &reader.functions {
        let offset = index_offset + output.len();
        output.extend(function_pseudo_ms2(function, offset)?);
    }
    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn groups_transitions_by_q1_in_first_seen_order() {
        let transitions = vec![
            TqMrmTransition {
                slot: 0,
                precursor_mz: 300.0,
                product_mz: 100.0,
            },
            TqMrmTransition {
                slot: 1,
                precursor_mz: 300.0,
                product_mz: 150.0,
            },
            TqMrmTransition {
                slot: 2,
                precursor_mz: 450.0,
                product_mz: 200.0,
            },
            TqMrmTransition {
                slot: 3,
                precursor_mz: 300.0,
                product_mz: 175.0,
            },
        ];

        let groups = transition_groups(&transitions);
        assert_eq!(groups.len(), 2);
        assert_eq!(groups[0].precursor_mz, 300.0);
        assert_eq!(groups[0].transition_indices, vec![0, 1, 3]);
        assert_eq!(groups[1].precursor_mz, 450.0);
        assert_eq!(groups[1].transition_indices, vec![2]);
    }

    #[test]
    fn summarize_sparse_transition_spectrum() {
        let (tic, bp_mz, bp_int, low, high) =
            summarize(&[100.0, 150.0, 175.0], &[2.0, 9.0, 4.0]);
        assert_eq!(tic, 15.0);
        assert_eq!(bp_mz, Some(150.0));
        assert_eq!(bp_int, Some(9.0));
        assert_eq!(low, Some(100.0));
        assert_eq!(high, Some(175.0));
    }
}
