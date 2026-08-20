//! PSI-mzML semantic corrections for the experimental Waters TQ export path.
//!
//! `openmassspec-core` 1.5.0 provides the generic mzML serializer used by the
//! established OpenWRaw readers. The TQ path needs a few stricter semantics
//! that are not yet expressible through that public API (dynamic fileContent,
//! a source-bundle checksum, generic/unknown dissociation, corrected units,
//! and independent observed-vs-declared Q3 mass ranges).
//! Rather than changing legacy output, TQ conversion serializes plain mzML to
//! memory, applies deterministic PSI corrections, and (when requested) builds
//! a fresh indexed-mzML wrapper from the corrected bytes.

use std::collections::BTreeMap;
use std::fs;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};

use openmassspec_core as msc;

use crate::raw::tq_reader::TqReader;

/// Serialize a TQ source as PSI-corrected, non-indexed mzML.
pub fn write_tq_psi_mzml<S: msc::SpectrumSource + ?Sized, P: AsRef<Path>, W: Write>(
    source: &mut S,
    source_dir: P,
    out: &mut W,
) -> crate::Result<()> {
    let corrected = corrected_plain_mzml(source, source_dir.as_ref())?;
    out.write_all(corrected.as_bytes())?;
    Ok(())
}

/// Serialize a TQ source as PSI-corrected indexed mzML.
///
/// Offsets are calculated only after semantic corrections have been applied.
/// The indexed wrapper checksum follows the mzML indexed schema definition:
/// SHA-1 from the beginning of the file through the end of the opening
/// `<fileChecksum>` tag.
pub fn write_tq_psi_indexed_mzml<
    S: msc::SpectrumSource + ?Sized,
    P: AsRef<Path>,
    W: Write,
>(
    source: &mut S,
    source_dir: P,
    out: &mut W,
) -> crate::Result<()> {
    let corrected = corrected_plain_mzml(source, source_dir.as_ref())?;
    let indexed = build_indexed_mzml(&corrected)?;
    out.write_all(&indexed)?;
    Ok(())
}

fn corrected_plain_mzml<S: msc::SpectrumSource + ?Sized>(
    source: &mut S,
    source_dir: &Path,
) -> crate::Result<String> {
    let mut raw = Vec::new();
    msc::write_mzml(source, &mut raw).map_err(crate::Error::Io)?;
    let xml = String::from_utf8(raw)
        .map_err(|err| crate::Error::Parse(format!("TQ mzML writer emitted non-UTF-8 XML: {err}")))?;
    let source_sha1 = deterministic_bundle_sha1(source_dir)?;
    let scan_windows = q3_scan_windows(source_dir)?;
    apply_psi_corrections(xml, &source_sha1, &scan_windows)
}

fn apply_psi_corrections(
    mut xml: String,
    source_sha1: &str,
    scan_windows: &BTreeMap<String, (f64, f64)>,
) -> crate::Result<String> {
    // fileContent must describe what is actually present in spectrumList,
    // rather than unconditionally advertising both MS1 and MSn.
    let spectrum_body = xml
        .split_once("<spectrumList")
        .map(|(_, body)| body)
        .unwrap_or("");
    let has_ms1 = spectrum_body.contains("accession=\"MS:1000579\"");
    let has_msn = spectrum_body.contains("accession=\"MS:1000580\"");

    let old_file_content = concat!(
        "      <cvParam cvRef=\"MS\" accession=\"MS:1000579\" name=\"MS1 spectrum\" value=\"\"/>\n",
        "      <cvParam cvRef=\"MS\" accession=\"MS:1000580\" name=\"MSn spectrum\" value=\"\"/>\n"
    );
    let mut new_file_content = String::new();
    if has_ms1 {
        new_file_content.push_str(
            "      <cvParam cvRef=\"MS\" accession=\"MS:1000579\" name=\"MS1 spectrum\" value=\"\"/>\n",
        );
    }
    if has_msn {
        new_file_content.push_str(
            "      <cvParam cvRef=\"MS\" accession=\"MS:1000580\" name=\"MSn spectrum\" value=\"\"/>\n",
        );
    }
    if !xml.contains(old_file_content) {
        return Err(crate::Error::Parse(
            "TQ mzML PSI correction: expected fileContent template was not found".to_owned(),
        ));
    }
    xml = xml.replacen(old_file_content, &new_file_content, 1);

    // Waters RAW is a directory bundle. The checksum convention is made
    // explicit so the SHA-1 is reproducible and not confused with a vendor
    // definition for hashing a directory object.
    let source_close = "      </sourceFile>";
    let source_insert = format!(
        concat!(
            "        <cvParam cvRef=\"MS\" accession=\"MS:1000569\" name=\"SHA-1\" value=\"{}\"/>\n",
            "        <userParam name=\"OpenWRaw Waters RAW bundle checksum convention\" ",
            "value=\"SHA-1 over sorted relative-path length/path, file-size, and file bytes; mzML outputs excluded\"/>\n",
            "{}"
        ),
        source_sha1, source_close
    );
    if !xml.contains(source_close) {
        return Err(crate::Error::Parse(
            "TQ mzML PSI correction: sourceFile closing element not found".to_owned(),
        ));
    }
    xml = xml.replacen(source_close, &source_insert, 1);

    // An MS2 precursor activation group must carry the dissociation-method
    // parent term or one of its children. When TQ method metadata has not yet
    // identified a specific mechanism, use the generic parent rather than
    // inventing CID.
    xml = xml.replace(
        "            <activation>\n            </activation>",
        concat!(
            "            <activation>\n",
            "              <cvParam cvRef=\"MS\" accession=\"MS:1000044\" name=\"dissociation method\" value=\"\"/>\n",
            "            </activation>"
        ),
    );

    // openmassspec-core 1.5.0 hardcodes CID on SRM chromatograms. Until .EE
    // method metadata is exposed, downgrade that assertion to the generic
    // dissociation-method term.
    xml = xml.replace(
        "<cvParam cvRef=\"MS\" accession=\"MS:1000133\" name=\"collision-induced dissociation\" value=\"\"/>",
        "<cvParam cvRef=\"MS\" accession=\"MS:1000044\" name=\"dissociation method\" value=\"\"/>",
    );

    // PSI units for spectrum-level mass/intensity attributes.
    xml = xml.replace(
        "name=\"base peak m/z\" value=\"",
        "name=\"base peak m/z\" unitCvRef=\"MS\" unitAccession=\"MS:1000040\" unitName=\"m/z\" value=\"",
    );
    xml = xml.replace(
        "name=\"base peak intensity\" value=\"",
        "name=\"base peak intensity\" unitCvRef=\"MS\" unitAccession=\"MS:1000131\" unitName=\"number of detector counts\" value=\"",
    );
    xml = xml.replace(
        "name=\"lowest observed m/z\" value=\"",
        "name=\"lowest observed m/z\" unitCvRef=\"MS\" unitAccession=\"MS:1000040\" unitName=\"m/z\" value=\"",
    );
    xml = xml.replace(
        "name=\"highest observed m/z\" value=\"",
        "name=\"highest observed m/z\" unitCvRef=\"MS\" unitAccession=\"MS:1000040\" unitName=\"m/z\" value=\"",
    );

    // Intensity arrays (spectra and chromatograms) should identify their unit.
    xml = xml.replace(
        "<cvParam cvRef=\"MS\" accession=\"MS:1000515\" name=\"intensity array\" value=\"\"/>",
        "<cvParam cvRef=\"MS\" accession=\"MS:1000515\" name=\"intensity array\" value=\"\" unitCvRef=\"MS\" unitAccession=\"MS:1000131\" unitName=\"number of detector counts\"/>",
    );

    // The generic serializer uses lowest/highest *observed* m/z as the scan
    // window. Q3 descriptors carry the programmed acquisition bounds
    // separately, so repair only scanWindow while leaving observed terms
    // untouched.
    apply_scan_window_overrides(&mut xml, scan_windows)?;

    Ok(xml)
}

fn q3_scan_windows(dir: &Path) -> crate::Result<BTreeMap<String, (f64, f64)>> {
    let reader = TqReader::open(dir)?;
    let mut windows = BTreeMap::new();
    for function in &reader.q3_functions {
        let low = function.descriptor.mz_low as f64;
        let high = function.descriptor.mz_high as f64;
        if !low.is_finite() || !high.is_finite() || low <= 0.0 || high <= low {
            continue;
        }
        for scan_index in 0..function.scan_count() {
            windows.insert(
                format!(
                    "function={} process=0 scan={}",
                    function.index,
                    scan_index + 1
                ),
                (low, high),
            );
        }
    }
    Ok(windows)
}

fn apply_scan_window_overrides(
    xml: &mut String,
    scan_windows: &BTreeMap<String, (f64, f64)>,
) -> crate::Result<()> {
    for (native_id, &(low, high)) in scan_windows {
        let marker = format!("<spectrum id=\"{native_id}\"");
        let Some(start) = xml.find(&marker) else {
            // A source adapter is allowed to omit an undecodable spectrum.
            // Its static RAW descriptor may therefore have no emitted XML
            // counterpart; this is not an error.
            continue;
        };
        let relative_end = xml[start..].find("</spectrum>").ok_or_else(|| {
            crate::Error::Parse(format!(
                "TQ mzML PSI correction: unterminated spectrum {native_id}"
            ))
        })?;
        let end = start + relative_end + "</spectrum>".len();
        let mut spectrum = xml[start..end].to_owned();
        spectrum = replace_cv_value(&spectrum, "MS:1000501", low)?;
        spectrum = replace_cv_value(&spectrum, "MS:1000500", high)?;
        xml.replace_range(start..end, &spectrum);
    }
    Ok(())
}

fn replace_cv_value(text: &str, accession: &str, value: f64) -> crate::Result<String> {
    let accession_marker = format!("accession=\"{accession}\"");
    let accession_start = text.find(&accession_marker).ok_or_else(|| {
        crate::Error::Parse(format!(
            "TQ mzML PSI correction: spectrum lacks {accession} scan-window term"
        ))
    })?;
    let tail = &text[accession_start..];
    let relative_value = tail.find("value=\"").ok_or_else(|| {
        crate::Error::Parse(format!(
            "TQ mzML PSI correction: {accession} term lacks a value attribute"
        ))
    })?;
    let value_start = accession_start + relative_value + "value=\"".len();
    let value_end = text[value_start..]
        .find('"')
        .map(|offset| value_start + offset)
        .ok_or_else(|| {
            crate::Error::Parse(format!(
                "TQ mzML PSI correction: {accession} value is unterminated"
            ))
        })?;

    let mut corrected = text.to_owned();
    corrected.replace_range(value_start..value_end, &format!("{value:.6}"));
    Ok(corrected)
}

fn build_indexed_mzml(plain: &str) -> crate::Result<Vec<u8>> {
    let first_newline = plain.find('\n').ok_or_else(|| {
        crate::Error::Parse("TQ indexed mzML: missing XML declaration newline".to_owned())
    })?;
    let declaration = &plain[..=first_newline];
    let mzml = &plain[first_newline + 1..];

    let mut out = Vec::new();
    out.extend_from_slice(declaration.as_bytes());
    out.extend_from_slice(
        concat!(
            "<indexedmzML xmlns=\"http://psi.hupo.org/ms/mzml\"\n",
            "             xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n",
            "             xsi:schemaLocation=\"http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.2_idx.xsd\">\n"
        )
        .as_bytes(),
    );
    out.extend_from_slice(mzml.as_bytes());

    let spectrum_offsets = collect_element_offsets(&out, b"<spectrum ", "id")?;
    let chromatogram_offsets = collect_element_offsets(&out, b"<chromatogram ", "id")?;

    let index_list_offset = out.len();
    let index_count = 1 + usize::from(!chromatogram_offsets.is_empty());
    writeln!(&mut out, "<indexList count=\"{index_count}\">")?;
    writeln!(&mut out, "  <index name=\"spectrum\">")?;
    for (id, offset) in &spectrum_offsets {
        writeln!(&mut out, "    <offset idRef=\"{id}\">{offset}</offset>")?;
    }
    writeln!(&mut out, "  </index>")?;
    if !chromatogram_offsets.is_empty() {
        writeln!(&mut out, "  <index name=\"chromatogram\">")?;
        for (id, offset) in &chromatogram_offsets {
            writeln!(&mut out, "    <offset idRef=\"{id}\">{offset}</offset>")?;
        }
        writeln!(&mut out, "  </index>")?;
    }
    writeln!(&mut out, "</indexList>")?;
    writeln!(
        &mut out,
        "<indexListOffset>{index_list_offset}</indexListOffset>"
    )?;

    out.extend_from_slice(b"<fileChecksum>");
    let digest = sha1_bytes(&out);
    let hex = hex_digest(&digest);
    out.extend_from_slice(hex.as_bytes());
    out.extend_from_slice(b"</fileChecksum>\n</indexedmzML>\n");
    Ok(out)
}

fn collect_element_offsets(
    bytes: &[u8],
    needle: &[u8],
    attribute: &str,
) -> crate::Result<Vec<(String, usize)>> {
    let mut result = Vec::new();
    let mut cursor = 0_usize;
    while let Some(relative) = find_bytes(&bytes[cursor..], needle) {
        let offset = cursor + relative;
        let end = bytes[offset..]
            .iter()
            .position(|&byte| byte == b'>')
            .map(|pos| offset + pos)
            .ok_or_else(|| crate::Error::Parse("TQ indexed mzML: unterminated element".to_owned()))?;
        let tag = std::str::from_utf8(&bytes[offset..=end])
            .map_err(|err| crate::Error::Parse(format!("TQ indexed mzML: invalid UTF-8 tag: {err}")))?;
        let id = extract_attribute(tag, attribute).ok_or_else(|| {
            crate::Error::Parse(format!(
                "TQ indexed mzML: element at byte {offset} lacks {attribute} attribute"
            ))
        })?;
        result.push((id.to_owned(), offset));
        cursor = end + 1;
    }
    Ok(result)
}

fn extract_attribute<'a>(tag: &'a str, name: &str) -> Option<&'a str> {
    let prefix = format!("{name}=\"");
    let start = tag.find(&prefix)? + prefix.len();
    let tail = &tag[start..];
    let end = tail.find('"')?;
    Some(&tail[..end])
}

fn find_bytes(haystack: &[u8], needle: &[u8]) -> Option<usize> {
    if needle.is_empty() {
        return Some(0);
    }
    haystack.windows(needle.len()).position(|window| window == needle)
}

fn deterministic_bundle_sha1(dir: &Path) -> crate::Result<String> {
    let mut files = Vec::new();
    collect_source_files(dir, dir, &mut files)?;
    files.sort_by(|left, right| left.0.cmp(&right.0));

    let mut sha = Sha1::new();
    for (relative, path) in files {
        let rel = relative.as_bytes();
        sha.update(&(rel.len() as u64).to_le_bytes());
        sha.update(rel);

        let size = fs::metadata(&path)?.len();
        sha.update(&size.to_le_bytes());

        let mut file = fs::File::open(path)?;
        let mut buffer = [0_u8; 64 * 1024];
        loop {
            let read = file.read(&mut buffer)?;
            if read == 0 {
                break;
            }
            sha.update(&buffer[..read]);
        }
    }
    Ok(hex_digest(&sha.finalize()))
}

fn collect_source_files(
    root: &Path,
    current: &Path,
    out: &mut Vec<(String, PathBuf)>,
) -> crate::Result<()> {
    for entry in fs::read_dir(current)? {
        let entry = entry?;
        let path = entry.path();
        let file_type = entry.file_type()?;
        if file_type.is_dir() {
            collect_source_files(root, &path, out)?;
            continue;
        }
        if !file_type.is_file() {
            continue;
        }

        let relative_path = path.strip_prefix(root).map_err(|err| {
            crate::Error::Parse(format!("TQ source checksum: cannot relativize path: {err}"))
        })?;
        let relative = relative_path.to_string_lossy().replace('\\', "/");
        let lower = relative.to_ascii_lowercase();
        if lower.ends_with(".mzml") || lower.ends_with(".mzml.gz") {
            continue;
        }
        out.push((relative, path));
    }
    Ok(())
}

fn sha1_bytes(bytes: &[u8]) -> [u8; 20] {
    let mut sha = Sha1::new();
    sha.update(bytes);
    sha.finalize()
}

fn hex_digest(digest: &[u8; 20]) -> String {
    digest.iter().map(|byte| format!("{byte:02x}")).collect()
}

struct Sha1 {
    state: [u32; 5],
    count: u64,
    buffer: [u8; 64],
    buffer_len: usize,
}

impl Sha1 {
    fn new() -> Self {
        Self {
            state: [0x6745_2301, 0xefcd_ab89, 0x98ba_dcfe, 0x1032_5476, 0xc3d2_e1f0],
            count: 0,
            buffer: [0; 64],
            buffer_len: 0,
        }
    }

    fn update(&mut self, bytes: &[u8]) {
        let mut offset = 0;
        while offset < bytes.len() {
            let available = 64 - self.buffer_len;
            let take = available.min(bytes.len() - offset);
            self.buffer[self.buffer_len..self.buffer_len + take]
                .copy_from_slice(&bytes[offset..offset + take]);
            self.buffer_len += take;
            self.count += take as u64;
            offset += take;
            if self.buffer_len == 64 {
                self.compress();
                self.buffer_len = 0;
            }
        }
    }

    fn compress(&mut self) {
        let mut words = [0_u32; 80];
        for (index, word) in words.iter_mut().enumerate().take(16) {
            let base = index * 4;
            *word = u32::from_be_bytes([
                self.buffer[base],
                self.buffer[base + 1],
                self.buffer[base + 2],
                self.buffer[base + 3],
            ]);
        }
        for index in 16..80 {
            words[index] = (words[index - 3]
                ^ words[index - 8]
                ^ words[index - 14]
                ^ words[index - 16])
                .rotate_left(1);
        }

        let [mut a, mut b, mut c, mut d, mut e] = self.state;
        for (index, &word) in words.iter().enumerate() {
            let (f, k) = match index {
                0..=19 => ((b & c) | ((!b) & d), 0x5a82_7999),
                20..=39 => (b ^ c ^ d, 0x6ed9_eba1),
                40..=59 => ((b & c) | (b & d) | (c & d), 0x8f1b_bcdc),
                _ => (b ^ c ^ d, 0xca62_c1d6),
            };
            let temp = a
                .rotate_left(5)
                .wrapping_add(f)
                .wrapping_add(e)
                .wrapping_add(k)
                .wrapping_add(word);
            e = d;
            d = c;
            c = b.rotate_left(30);
            b = a;
            a = temp;
        }

        self.state[0] = self.state[0].wrapping_add(a);
        self.state[1] = self.state[1].wrapping_add(b);
        self.state[2] = self.state[2].wrapping_add(c);
        self.state[3] = self.state[3].wrapping_add(d);
        self.state[4] = self.state[4].wrapping_add(e);
    }

    fn finalize(mut self) -> [u8; 20] {
        let bit_count = self.count * 8;
        self.update(&[0x80]);
        while self.buffer_len != 56 {
            self.update(&[0]);
        }
        self.update(&bit_count.to_be_bytes());

        let mut digest = [0_u8; 20];
        for (index, word) in self.state.iter().enumerate() {
            digest[index * 4..index * 4 + 4].copy_from_slice(&word.to_be_bytes());
        }
        digest
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sha1_matches_known_vector() {
        assert_eq!(
            hex_digest(&sha1_bytes(b"abc")),
            "a9993e364706816aba3e25717850c26c9cd0d89d"
        );
    }

    #[test]
    fn corrections_make_file_content_dynamic_and_add_units() {
        let xml = concat!(
            "<fileContent>\n",
            "      <cvParam cvRef=\"MS\" accession=\"MS:1000579\" name=\"MS1 spectrum\" value=\"\"/>\n",
            "      <cvParam cvRef=\"MS\" accession=\"MS:1000580\" name=\"MSn spectrum\" value=\"\"/>\n",
            "</fileContent>\n",
            "<sourceFile>\n      </sourceFile>\n",
            "<spectrumList><spectrum><cvParam cvRef=\"MS\" accession=\"MS:1000579\" name=\"MS1 spectrum\" value=\"\"/>",
            "<cvParam cvRef=\"MS\" accession=\"MS:1000504\" name=\"base peak m/z\" value=\"100\"/>",
            "<binaryDataArray><cvParam cvRef=\"MS\" accession=\"MS:1000515\" name=\"intensity array\" value=\"\"/></binaryDataArray>",
            "</spectrum></spectrumList>"
        );
        let fixed = apply_psi_corrections(xml.to_owned(), "deadbeef", &BTreeMap::new()).unwrap();
        let header = fixed.split("<spectrumList").next().unwrap();
        assert!(header.contains("MS:1000579"));
        assert!(!header.contains("MS:1000580"));
        assert!(fixed.contains("unitAccession=\"MS:1000040\""));
        assert!(fixed.contains("unitAccession=\"MS:1000131\""));
        assert!(fixed.contains("MS:1000569"));
    }

    #[test]
    fn declared_scan_window_replaces_only_scan_window_values() {
        let mut xml = concat!(
            "<spectrum id=\"function=2 process=0 scan=1\">",
            "<cvParam accession=\"MS:1000528\" name=\"lowest observed m/z\" value=\"100.000000\"/>",
            "<cvParam accession=\"MS:1000527\" name=\"highest observed m/z\" value=\"200.000000\"/>",
            "<scanWindow>",
            "<cvParam accession=\"MS:1000501\" name=\"scan window lower limit\" value=\"100.000000\"/>",
            "<cvParam accession=\"MS:1000500\" name=\"scan window upper limit\" value=\"200.000000\"/>",
            "</scanWindow></spectrum>"
        )
        .to_owned();
        let mut windows = BTreeMap::new();
        windows.insert(
            "function=2 process=0 scan=1".to_owned(),
            (75.0, 900.0),
        );
        apply_scan_window_overrides(&mut xml, &windows).unwrap();
        assert!(xml.contains("lowest observed m/z\" value=\"100.000000\""));
        assert!(xml.contains("highest observed m/z\" value=\"200.000000\""));
        assert!(xml.contains("scan window lower limit\" value=\"75.000000\""));
        assert!(xml.contains("scan window upper limit\" value=\"900.000000\""));
    }

    #[test]
    fn indexed_checksum_includes_open_file_checksum_tag() {
        let plain = concat!(
            "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n",
            "<mzML xmlns=\"http://psi.hupo.org/ms/mzml\" version=\"1.1.0\">\n",
            "<run><spectrumList count=\"1\"><spectrum id=\"function=1 process=0 scan=1\" index=\"0\" defaultArrayLength=\"0\"></spectrum></spectrumList></run>\n",
            "</mzML>\n"
        );
        let indexed = build_indexed_mzml(plain).unwrap();
        let text = String::from_utf8(indexed.clone()).unwrap();
        let index_offset_start = text.find("<indexListOffset>").unwrap() + "<indexListOffset>".len();
        let index_offset_end = text[index_offset_start..].find('<').unwrap() + index_offset_start;
        let index_offset: usize = text[index_offset_start..index_offset_end].parse().unwrap();
        assert!(text[index_offset..].starts_with("<indexList"));

        let checksum_start = text.find("<fileChecksum>").unwrap();
        let value_start = checksum_start + "<fileChecksum>".len();
        let value_end = text[value_start..].find('<').unwrap() + value_start;
        let expected = hex_digest(&sha1_bytes(&indexed[..value_start]));
        assert_eq!(&text[value_start..value_end], expected);
    }
}
