use openmassspec_core::SpectrumSource;
use openwraw::raw::tq::{TqPolarity, FUNCTION_RECORD_SIZE};
use openwraw::raw::tq_mrm::TqMrmReader;
use openwraw::raw::tq_reader::TqReader;
use openwraw::tq_mixed_mzml::{write_tq_mixed_mzml, TqMixedSource};
use openwraw::tq_mzml::TqQ3MzmlMode;
use std::fs;
use std::path::PathBuf;

fn temp_bundle(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("openwraw-{name}-{}", std::process::id()))
}

fn idx_record(offset: u32, format_marker: u32, count: u32, rt_min: f32) -> [u8; 22] {
    let mut record = [0_u8; 22];
    record[0..4].copy_from_slice(&offset.to_le_bytes());
    record[4..8].copy_from_slice(&(format_marker | count).to_le_bytes());
    record[12..16].copy_from_slice(&rt_min.to_le_bytes());
    record
}

fn packed_mrm_value(power: u32, base: u32) -> [u8; 4] {
    ((power << 22) | (base & 0x001f_ffff)).to_le_bytes()
}

fn direct6_record(
    mass_base: u32,
    mass_power_field: u32,
    intensity_base: i16,
    intensity_power: u32,
) -> [u8; 6] {
    let packed = (mass_base << 9)
        | ((mass_power_field & 0x1f) << 4)
        | (intensity_power & 0x0f);
    let mut record = [0_u8; 6];
    record[0..2].copy_from_slice(&intensity_base.to_le_bytes());
    record[2..6].copy_from_slice(&packed.to_le_bytes());
    record
}

fn write_synthetic_mixed_bundle(dir: &PathBuf) {
    let _ = fs::remove_dir_all(dir);
    fs::create_dir_all(dir).unwrap();

    fs::write(
        dir.join("_HEADER.TXT"),
        "$$ Version: 01.00\r\n\
$$ Instrument: SYNTHETIC-TQ\r\n\
$$ Cal Function 2: 0.0,1.0,T0\r\n",
    )
    .unwrap();

    let mut functions = vec![0_u8; FUNCTION_RECORD_SIZE * 2];

    // Function 1: positive MRM with two synthetic transitions.
    functions[0] = 0x09;
    functions[0x0a0..0x0a4].copy_from_slice(&300.0_f32.to_le_bytes());
    functions[0x0a4..0x0a8].copy_from_slice(&300.0_f32.to_le_bytes());
    functions[0x120..0x124].copy_from_slice(&100.0_f32.to_le_bytes());
    functions[0x124..0x128].copy_from_slice(&150.0_f32.to_le_bytes());

    // Function 2: negative broad Q3 scan.
    let f2 = FUNCTION_RECORD_SIZE;
    functions[f2] = 0x2b;
    functions[f2 + 0x018..f2 + 0x01c].copy_from_slice(&40.0_f32.to_le_bytes());
    functions[f2 + 0x020..f2 + 0x024].copy_from_slice(&0.5_f32.to_le_bytes());
    functions[f2 + 0x0a0..f2 + 0x0a4].copy_from_slice(&75.0_f32.to_le_bytes());
    functions[f2 + 0x120..f2 + 0x124].copy_from_slice(&900.0_f32.to_le_bytes());
    fs::write(dir.join("_FUNCTNS.INF"), functions).unwrap();

    // MRM: two cycles, two channels per cycle.
    let mut mrm_idx = Vec::new();
    mrm_idx.extend_from_slice(&idx_record(0, 0x0800_0000, 2, 0.0));
    mrm_idx.extend_from_slice(&idx_record(8, 0x0800_0000, 2, 0.5));
    fs::write(dir.join("_FUNC001.IDX"), mrm_idx).unwrap();

    let mut mrm_dat = Vec::new();
    mrm_dat.extend_from_slice(&packed_mrm_value(10, 1024)); // 1
    mrm_dat.extend_from_slice(&packed_mrm_value(11, 1024)); // 2
    mrm_dat.extend_from_slice(&packed_mrm_value(12, 1024)); // 4
    mrm_dat.extend_from_slice(&packed_mrm_value(13, 1024)); // 8
    fs::write(dir.join("_FUNC001.DAT"), mrm_dat).unwrap();

    // Q3: one direct-6 profile scan with two points.
    let q3_idx = idx_record(0, 0x1800_0000, 2, 0.25);
    fs::write(dir.join("_FUNC002.IDX"), q3_idx).unwrap();
    let mut q3_dat = Vec::new();
    q3_dat.extend_from_slice(&direct6_record(100, 23, 10, 0));
    q3_dat.extend_from_slice(&direct6_record(201, 22, 5, 1));
    fs::write(dir.join("_FUNC002.DAT"), q3_dat).unwrap();
}

#[test]
fn synthetic_mixed_bundle_decodes_q3_and_mrm_independently() {
    let dir = temp_bundle("tq-mixed-decode");
    write_synthetic_mixed_bundle(&dir);

    let q3 = TqReader::open(&dir).unwrap();
    assert_eq!(q3.q3_functions.len(), 1);
    assert_eq!(q3.mrm_function_count, 1);
    let scan = q3.decode_scan(2, 0).unwrap();
    assert_eq!(scan.polarity, Some(TqPolarity::Negative));
    assert_eq!(scan.spectrum.mz, vec![100.0, 100.5]);
    assert_eq!(scan.spectrum.intensity, vec![10.0, 20.0]);

    let mrm = TqMrmReader::open(&dir).unwrap();
    assert_eq!(mrm.functions.len(), 1);
    assert_eq!(mrm.transition_count(), 2);
    let traces = mrm.chromatograms().unwrap();
    assert_eq!(traces.len(), 2);
    assert_eq!(traces[0].precursor_mz, 300.0);
    assert_eq!(traces[0].product_mz, 100.0);
    assert_eq!(traces[0].time_sec, vec![0.0, 30.0]);
    assert_eq!(traces[0].intensity, vec![1.0, 4.0]);
    assert_eq!(traces[1].intensity, vec![2.0, 8.0]);

    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn mixed_source_emits_q3_spectrum_and_srm_chromatograms() {
    let dir = temp_bundle("tq-mixed-source");
    write_synthetic_mixed_bundle(&dir);

    let mut source = TqMixedSource::open(&dir, TqQ3MzmlMode::PseudoMs1).unwrap();
    let spectra: Vec<_> = source.iter_spectra().collect();
    assert_eq!(spectra.len(), 1);
    assert_eq!(spectra[0].ms_level, 1);
    assert_eq!(spectra[0].polarity, Some(openmassspec_core::Polarity::Negative));

    let chromatograms: Vec<_> = source.iter_chromatograms().collect();
    assert_eq!(chromatograms.len(), 2);
    for chromatogram in &chromatograms {
        let cv = chromatogram.chromatogram_type.as_ref().unwrap();
        assert_eq!(cv.accession, "MS:1001473");
        assert!(chromatogram.precursor_mz.is_some());
        assert!(chromatogram.product_mz.is_some());
    }

    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn mixed_mzml_contains_spectrum_and_chromatogram_lists() {
    let dir = temp_bundle("tq-mixed-mzml");
    write_synthetic_mixed_bundle(&dir);

    let mut output = Vec::new();
    write_tq_mixed_mzml(&dir, &mut output, TqQ3MzmlMode::PseudoMs1).unwrap();
    let xml = String::from_utf8(output).unwrap();

    assert!(xml.contains("<spectrumList"));
    assert!(xml.contains("<chromatogramList"));
    assert!(xml.contains("selected reaction monitoring chromatogram"));

    let _ = fs::remove_dir_all(&dir);
}
