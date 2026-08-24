use openmassspec_core::SpectrumSource;
use openwraw::raw::tq::FUNCTION_RECORD_SIZE;
use openwraw::tq_mzml::{TqQ3MzmlMode, TqQ3Source};
use std::fs;
use std::path::PathBuf;

fn temp_bundle() -> PathBuf {
    std::env::temp_dir().join(format!("openwraw-tq-q3-failfast-{}", std::process::id()))
}

fn idx_record(offset: u32, count: u32, rt_min: f32) -> [u8; 22] {
    let mut record = [0_u8; 22];
    record[0..4].copy_from_slice(&offset.to_le_bytes());
    record[4..8].copy_from_slice(&(0x1800_0000_u32 | count).to_le_bytes());
    record[12..16].copy_from_slice(&rt_min.to_le_bytes());
    record
}

fn direct6_record() -> [u8; 6] {
    let mass_base = 100_u32;
    let packed = (mass_base << 9) | (23_u32 << 4);
    let mut record = [0_u8; 6];
    record[0..2].copy_from_slice(&10_i16.to_le_bytes());
    record[2..6].copy_from_slice(&packed.to_le_bytes());
    record
}

fn write_bundle(dir: &PathBuf) {
    let _ = fs::remove_dir_all(dir);
    fs::create_dir_all(dir).unwrap();
    fs::write(
        dir.join("_HEADER.TXT"),
        "$$ Version: 01.00\r\n$$ Instrument: SYNTHETIC-TQ\r\n$$ Cal Function 1: 0.0,1.0,T0\r\n",
    )
    .unwrap();

    let mut function = vec![0_u8; FUNCTION_RECORD_SIZE];
    function[0] = 0x0b;
    function[0x018..0x01c].copy_from_slice(&50.0_f32.to_le_bytes());
    function[0x020..0x024].copy_from_slice(&1.0_f32.to_le_bytes());
    function[0x0a0..0x0a4].copy_from_slice(&90.0_f32.to_le_bytes());
    function[0x120..0x124].copy_from_slice(&1350.0_f32.to_le_bytes());
    fs::write(dir.join("_FUNCTNS.INF"), function).unwrap();
    fs::write(dir.join("_FUNC001.IDX"), idx_record(0, 1, 0.1)).unwrap();
    fs::write(dir.join("_FUNC001.DAT"), direct6_record()).unwrap();
}

#[test]
fn runtime_q3_io_failure_is_fail_fast_not_silent() {
    let dir = temp_bundle();
    write_bundle(&dir);

    let mut source = TqQ3Source::open(&dir, TqQ3MzmlMode::PseudoMs1).unwrap();
    fs::remove_file(dir.join("_FUNC001.DAT")).unwrap();

    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        source.iter_spectra().collect::<Vec<_>>()
    }));
    assert!(
        result.is_err(),
        "runtime Q3 I/O failures must not silently reduce the spectrum stream"
    );

    let _ = fs::remove_dir_all(&dir);
}
