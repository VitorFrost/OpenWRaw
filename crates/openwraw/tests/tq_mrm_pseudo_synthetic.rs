use openmassspec_core::conformance::assert_source_invariants;
use openwraw::raw::tq::FUNCTION_RECORD_SIZE;
use openwraw::raw::tq_mrm::TqMrmReader;
use openwraw::tq_mixed_mzml::TqMixedSource;
use openwraw::tq_mrm_spectra::{pseudo_ms2_records, TqMrmSpectrumMode};
use openwraw::tq_mzml::TqQ3MzmlMode;
use std::fs;
use std::path::{Path, PathBuf};

fn temp_bundle() -> PathBuf {
    std::env::temp_dir().join(format!("openwraw-tq-mrm-pseudo-{}", std::process::id()))
}

fn idx_record(offset: u32, count: u32, rt_min: f32) -> [u8; 22] {
    let mut record = [0_u8; 22];
    record[0..4].copy_from_slice(&offset.to_le_bytes());
    record[4..8].copy_from_slice(&(0x0800_0000_u32 | count).to_le_bytes());
    record[12..16].copy_from_slice(&rt_min.to_le_bytes());
    record
}

fn packed_value(power: u32, base: u32) -> [u8; 4] {
    ((power << 22) | (base & 0x001f_ffff)).to_le_bytes()
}

fn write_bundle(dir: &Path) {
    let _ = fs::remove_dir_all(dir);
    fs::create_dir_all(dir).unwrap();
    fs::write(dir.join("_HEADER.TXT"), "$$ Version: 01.00\r\n").unwrap();

    // One positive MRM function with 3 transitions split across 2 synthetic Q1 groups.
    // The first Q1 group's Q3 values are deliberately stored out of m/z order
    // to verify the pseudo-spectrum projection sorts (Q3, intensity) pairs together.
    let mut function = [0_u8; FUNCTION_RECORD_SIZE];
    function[0] = 0x09;
    function[0x0a0..0x0a4].copy_from_slice(&300.0_f32.to_le_bytes());
    function[0x0a4..0x0a8].copy_from_slice(&300.0_f32.to_le_bytes());
    function[0x0a8..0x0ac].copy_from_slice(&450.0_f32.to_le_bytes());
    function[0x120..0x124].copy_from_slice(&150.0_f32.to_le_bytes());
    function[0x124..0x128].copy_from_slice(&100.0_f32.to_le_bytes());
    function[0x128..0x12c].copy_from_slice(&200.0_f32.to_le_bytes());
    fs::write(dir.join("_FUNCTNS.INF"), function).unwrap();

    let mut idx = Vec::new();
    idx.extend_from_slice(&idx_record(0, 3, 0.25));
    idx.extend_from_slice(&idx_record(12, 3, 0.50));
    fs::write(dir.join("_FUNC001.IDX"), idx).unwrap();

    let mut dat = Vec::new();
    // Cycle 1, transition order: Q3=150 -> 1; Q3=100 -> 2; Q3=200 -> 4.
    dat.extend_from_slice(&packed_value(10, 1024));
    dat.extend_from_slice(&packed_value(11, 1024));
    dat.extend_from_slice(&packed_value(12, 1024));
    // Cycle 2: Q3=150 -> 8; Q3=100 -> 16; Q3=200 -> 32.
    dat.extend_from_slice(&packed_value(13, 1024));
    dat.extend_from_slice(&packed_value(14, 1024));
    dat.extend_from_slice(&packed_value(15, 1024));
    fs::write(dir.join("_FUNC001.DAT"), dat).unwrap();
}

#[test]
fn pseudo_ms2_groups_by_q1_sorts_q3_and_preserves_signal_pairing() {
    let dir = temp_bundle();
    write_bundle(&dir);

    let reader = TqMrmReader::open(&dir).unwrap();
    let spectra = pseudo_ms2_records(&reader, 7).unwrap();

    // 2 cycles × 2 Q1 groups = 4 sparse pseudo-MS2 spectra.
    assert_eq!(spectra.len(), 4);
    assert_eq!(spectra[0].index, 7);
    assert_eq!(spectra[3].index, 10);
    assert_eq!(spectra[0].scan_number, 8);
    assert_eq!(spectra[0].ms_level, 2);

    let p0 = spectra[0].precursor.as_ref().unwrap();
    assert_eq!(p0.target_mz, Some(300.0));
    // Q3 is sorted ascending and the signal values move with their Q3 values.
    assert_eq!(spectra[0].mz, vec![100.0, 150.0]);
    assert_eq!(spectra[0].intensity, vec![2.0, 1.0]);
    assert_eq!(spectra[0].retention_time_sec, 15.0);

    let p1 = spectra[1].precursor.as_ref().unwrap();
    assert_eq!(p1.target_mz, Some(450.0));
    assert_eq!(spectra[1].mz, vec![200.0]);
    assert_eq!(spectra[1].intensity, vec![4.0]);

    // Second cycle retains the same sorted Q1/Q3 grouping with new signal values.
    assert_eq!(spectra[2].mz, vec![100.0, 150.0]);
    assert_eq!(spectra[2].intensity, vec![16.0, 8.0]);
    assert_eq!(spectra[3].intensity, vec![32.0]);
    assert_eq!(spectra[2].retention_time_sec, 30.0);

    assert!(spectra.iter().all(|s| {
        s.filter
            .as_deref()
            .unwrap_or_default()
            .contains("pseudo-MS2")
    }));

    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn pseudo_ms2_mixed_source_satisfies_openmassspec_conformance_contract() {
    let dir = temp_bundle();
    write_bundle(&dir);

    // This synthetic bundle intentionally has no Q3 function, so the mixed
    // source contains exactly the four MRM-derived pseudo-MS2 spectra above.
    let mut source = TqMixedSource::open_with_mrm_spectra(
        &dir,
        TqQ3MzmlMode::NativeMs2,
        TqMrmSpectrumMode::PseudoMs2,
    )
    .unwrap();

    let count = assert_source_invariants(&mut source).unwrap();
    assert_eq!(count, 4);

    let _ = fs::remove_dir_all(&dir);
}
