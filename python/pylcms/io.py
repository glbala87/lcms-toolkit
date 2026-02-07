"""
I/O functions for reading and writing LC-MS data files.
"""

import base64
import struct
import zlib
from typing import Optional, List, Dict, Any, Callable
from pathlib import Path
import xml.etree.ElementTree as ET

from .spectrum import Spectrum, SpectrumType, Polarity, Precursor
from .chromatogram import Chromatogram, ChromatogramType
from .experiment import MSExperiment


# CV term accessions
CV_MS_LEVEL = "MS:1000511"
CV_PROFILE = "MS:1000128"
CV_CENTROID = "MS:1000127"
CV_POSITIVE = "MS:1000130"
CV_NEGATIVE = "MS:1000129"
CV_MZ_ARRAY = "MS:1000514"
CV_INTENSITY_ARRAY = "MS:1000515"
CV_TIME_ARRAY = "MS:1000595"
CV_FLOAT64 = "MS:1000523"
CV_FLOAT32 = "MS:1000521"
CV_ZLIB = "MS:1000574"
CV_SCAN_TIME = "MS:1000016"
CV_TIC = "MS:1000285"
CV_SELECTED_MZ = "MS:1000744"
CV_CHARGE = "MS:1000041"
CV_COLLISION_ENERGY = "MS:1000045"
CV_MINUTE = "UO:0000031"


def decode_binary(
    data: str,
    is_64bit: bool = True,
    is_compressed: bool = False,
    is_little_endian: bool = True,
) -> List[float]:
    """
    Decode Base64 encoded binary data.

    Args:
        data: Base64 encoded string
        is_64bit: True for 64-bit floats, False for 32-bit
        is_compressed: True if zlib compressed
        is_little_endian: True for little-endian byte order

    Returns:
        List of decoded values
    """
    # Remove whitespace
    data = "".join(data.split())
    if not data:
        return []

    # Decode Base64
    binary = base64.b64decode(data)

    # Decompress if needed
    if is_compressed:
        binary = zlib.decompress(binary)

    # Unpack floats
    if is_64bit:
        fmt = "<d" if is_little_endian else ">d"
        size = 8
    else:
        fmt = "<f" if is_little_endian else ">f"
        size = 4

    count = len(binary) // size
    values = []
    for i in range(count):
        value = struct.unpack(fmt, binary[i * size : (i + 1) * size])[0]
        values.append(value)

    return values


def load_mzml(
    filename: str,
    ms_levels: Optional[List[int]] = None,
    rt_range: Optional[tuple] = None,
    max_spectra: int = 0,
    skip_chromatograms: bool = False,
    progress_callback: Optional[Callable[[int, int], bool]] = None,
) -> MSExperiment:
    """
    Load an mzML file.

    Args:
        filename: Path to mzML file
        ms_levels: Only load spectra at these MS levels (None = all)
        rt_range: Only load spectra in (min, max) RT range
        max_spectra: Maximum number of spectra to load (0 = unlimited)
        skip_chromatograms: Skip loading chromatograms
        progress_callback: Callback(current, total) -> continue

    Returns:
        Loaded MSExperiment

    Example:
        >>> exp = load_mzml("sample.mzML")
        >>> exp = load_mzml("sample.mzML", ms_levels=[1], rt_range=(60, 300))
    """
    exp = MSExperiment()
    exp.source_file = filename

    # Parse XML
    tree = ET.parse(filename)
    root = tree.getroot()

    # Handle namespace
    ns = {}
    if root.tag.startswith("{"):
        ns_uri = root.tag.split("}")[0][1:]
        ns = {"mzml": ns_uri}

    # Find mzML element (may be inside indexedmzML)
    mzml = root.find(".//mzml:mzML", ns) or root.find(".//mzML", ns) or root

    # Find run element
    run = (
        mzml.find("mzml:run", ns)
        or mzml.find("run", ns)
        or mzml.find(".//run")
    )
    if run is None:
        run = mzml

    # Parse spectra
    spectrum_list = (
        run.find("mzml:spectrumList", ns)
        or run.find("spectrumList", ns)
        or run.find(".//spectrumList")
    )

    if spectrum_list is not None:
        total = int(spectrum_list.get("count", 0))
        loaded = 0

        for spec_elem in spectrum_list.findall("mzml:spectrum", ns) or spectrum_list.findall("spectrum"):
            if max_spectra > 0 and loaded >= max_spectra:
                break

            spec = _parse_spectrum(spec_elem, ns)

            # Apply filters
            if ms_levels is not None and spec.ms_level not in ms_levels:
                continue
            if rt_range is not None:
                if spec.rt < rt_range[0] or spec.rt > rt_range[1]:
                    continue

            exp.add_spectrum(spec)
            loaded += 1

            if progress_callback:
                if not progress_callback(loaded, total):
                    break

    # Parse chromatograms
    if not skip_chromatograms:
        chrom_list = (
            run.find("mzml:chromatogramList", ns)
            or run.find("chromatogramList", ns)
            or run.find(".//chromatogramList")
        )

        if chrom_list is not None:
            for chrom_elem in chrom_list.findall("mzml:chromatogram", ns) or chrom_list.findall("chromatogram"):
                chrom = _parse_chromatogram(chrom_elem, ns)
                exp.add_chromatogram(chrom)

    return exp


def _parse_spectrum(elem: ET.Element, ns: Dict[str, str]) -> Spectrum:
    """Parse a spectrum element."""
    spec = Spectrum()
    spec.native_id = elem.get("id", "")

    # Parse cvParams
    for cv in elem.findall("mzml:cvParam", ns) or elem.findall("cvParam"):
        acc = cv.get("accession", "")
        value = cv.get("value", "")

        if acc == CV_MS_LEVEL:
            spec.ms_level = int(value)
        elif acc == CV_PROFILE:
            spec.spectrum_type = SpectrumType.PROFILE
        elif acc == CV_CENTROID:
            spec.spectrum_type = SpectrumType.CENTROID
        elif acc == CV_POSITIVE:
            spec.polarity = Polarity.POSITIVE
        elif acc == CV_NEGATIVE:
            spec.polarity = Polarity.NEGATIVE

    # Parse scan list for RT
    scan_list = elem.find("mzml:scanList", ns) or elem.find("scanList")
    if scan_list is not None:
        scan = scan_list.find("mzml:scan", ns) or scan_list.find("scan")
        if scan is not None:
            for cv in scan.findall("mzml:cvParam", ns) or scan.findall("cvParam"):
                if cv.get("accession") == CV_SCAN_TIME:
                    rt = float(cv.get("value", 0))
                    if cv.get("unitAccession") == CV_MINUTE:
                        rt *= 60.0
                    spec.rt = rt

    # Parse precursors
    prec_list = elem.find("mzml:precursorList", ns) or elem.find("precursorList")
    if prec_list is not None:
        for prec_elem in prec_list.findall("mzml:precursor", ns) or prec_list.findall("precursor"):
            prec = _parse_precursor(prec_elem, ns)
            spec.precursors.append(prec)

    # Parse binary data
    binary_list = elem.find("mzml:binaryDataArrayList", ns) or elem.find("binaryDataArrayList")
    if binary_list is not None:
        mz_data = []
        intensity_data = []

        for binary in binary_list.findall("mzml:binaryDataArray", ns) or binary_list.findall("binaryDataArray"):
            is_mz = False
            is_intensity = False
            is_64bit = True
            is_compressed = False

            for cv in binary.findall("mzml:cvParam", ns) or binary.findall("cvParam"):
                acc = cv.get("accession", "")
                if acc == CV_MZ_ARRAY:
                    is_mz = True
                elif acc == CV_INTENSITY_ARRAY:
                    is_intensity = True
                elif acc == CV_FLOAT64:
                    is_64bit = True
                elif acc == CV_FLOAT32:
                    is_64bit = False
                elif acc == CV_ZLIB:
                    is_compressed = True

            data_elem = binary.find("mzml:binary", ns) or binary.find("binary")
            if data_elem is not None and data_elem.text:
                values = decode_binary(data_elem.text, is_64bit, is_compressed)
                if is_mz:
                    mz_data = values
                elif is_intensity:
                    intensity_data = values

        if mz_data:
            if not intensity_data:
                intensity_data = [0.0] * len(mz_data)
            spec.mz = mz_data
            spec.intensity = intensity_data

    return spec


def _parse_precursor(elem: ET.Element, ns: Dict[str, str]) -> Precursor:
    """Parse a precursor element."""
    prec = Precursor()

    # Parse isolation window
    isolation = elem.find("mzml:isolationWindow", ns) or elem.find("isolationWindow")
    if isolation is not None:
        for cv in isolation.findall("mzml:cvParam", ns) or isolation.findall("cvParam"):
            acc = cv.get("accession", "")
            value = float(cv.get("value", 0))
            if "1000827" in acc:  # target
                prec.mz = value
            elif "1000828" in acc:  # lower offset
                prec.isolation_window_lower = value
            elif "1000829" in acc:  # upper offset
                prec.isolation_window_upper = value

    # Parse selected ion
    selected_list = elem.find("mzml:selectedIonList", ns) or elem.find("selectedIonList")
    if selected_list is not None:
        selected = selected_list.find("mzml:selectedIon", ns) or selected_list.find("selectedIon")
        if selected is not None:
            for cv in selected.findall("mzml:cvParam", ns) or selected.findall("cvParam"):
                acc = cv.get("accession", "")
                value = cv.get("value", "0")
                if acc == CV_SELECTED_MZ:
                    prec.mz = float(value)
                elif acc == CV_CHARGE:
                    prec.charge = int(value)
                elif "1000042" in acc:  # intensity
                    prec.intensity = float(value)

    # Parse activation
    activation = elem.find("mzml:activation", ns) or elem.find("activation")
    if activation is not None:
        for cv in activation.findall("mzml:cvParam", ns) or activation.findall("cvParam"):
            acc = cv.get("accession", "")
            if acc == CV_COLLISION_ENERGY:
                prec.collision_energy = float(cv.get("value", 0))
            elif "1000133" in acc:
                prec.activation_method = "CID"
            elif "1000422" in acc:
                prec.activation_method = "HCD"
            elif "1000598" in acc:
                prec.activation_method = "ETD"

    return prec


def _parse_chromatogram(elem: ET.Element, ns: Dict[str, str]) -> Chromatogram:
    """Parse a chromatogram element."""
    chrom = Chromatogram()
    chrom.native_id = elem.get("id", "")

    # Parse cvParams
    for cv in elem.findall("mzml:cvParam", ns) or elem.findall("cvParam"):
        acc = cv.get("accession", "")
        if "1000235" in acc:  # TIC
            chrom.chrom_type = ChromatogramType.TIC
        elif "1000628" in acc:  # BPC
            chrom.chrom_type = ChromatogramType.BPC
        elif "1001473" in acc:  # SRM
            chrom.chrom_type = ChromatogramType.SRM

    # Parse binary data
    binary_list = elem.find("mzml:binaryDataArrayList", ns) or elem.find("binaryDataArrayList")
    if binary_list is not None:
        rt_data = []
        intensity_data = []

        for binary in binary_list.findall("mzml:binaryDataArray", ns) or binary_list.findall("binaryDataArray"):
            is_time = False
            is_intensity = False
            is_64bit = True
            is_compressed = False
            is_minutes = False

            for cv in binary.findall("mzml:cvParam", ns) or binary.findall("cvParam"):
                acc = cv.get("accession", "")
                unit = cv.get("unitAccession", "")
                if acc == CV_TIME_ARRAY:
                    is_time = True
                elif acc == CV_INTENSITY_ARRAY:
                    is_intensity = True
                elif acc == CV_FLOAT64:
                    is_64bit = True
                elif acc == CV_FLOAT32:
                    is_64bit = False
                elif acc == CV_ZLIB:
                    is_compressed = True
                if unit == CV_MINUTE:
                    is_minutes = True

            data_elem = binary.find("mzml:binary", ns) or binary.find("binary")
            if data_elem is not None and data_elem.text:
                values = decode_binary(data_elem.text, is_64bit, is_compressed)
                if is_time:
                    if is_minutes:
                        values = [v * 60.0 for v in values]
                    rt_data = values
                elif is_intensity:
                    intensity_data = values

        if rt_data and intensity_data:
            chrom.rt = rt_data
            chrom.intensity = intensity_data

    return chrom


def load_mzxml(
    filename: str,
    ms_levels: Optional[List[int]] = None,
    rt_range: Optional[tuple] = None,
    max_spectra: int = 0,
    progress_callback: Optional[Callable[[int, int], bool]] = None,
) -> MSExperiment:
    """
    Load an mzXML file.

    Args:
        filename: Path to mzXML file
        ms_levels: Only load spectra at these MS levels
        rt_range: Only load spectra in (min, max) RT range
        max_spectra: Maximum number of spectra to load
        progress_callback: Callback for progress updates

    Returns:
        Loaded MSExperiment
    """
    exp = MSExperiment()
    exp.source_file = filename

    tree = ET.parse(filename)
    root = tree.getroot()

    # Handle namespace
    ns = {}
    if root.tag.startswith("{"):
        ns_uri = root.tag.split("}")[0][1:]
        ns = {"mzxml": ns_uri}

    # Find msRun
    ms_run = root.find(".//mzxml:msRun", ns) or root.find(".//msRun") or root

    loaded = 0

    def parse_scans(parent):
        nonlocal loaded
        for scan in parent.findall("mzxml:scan", ns) or parent.findall("scan"):
            if max_spectra > 0 and loaded >= max_spectra:
                return

            spec = _parse_mzxml_scan(scan, ns)

            # Apply filters
            if ms_levels is not None and spec.ms_level not in ms_levels:
                parse_scans(scan)
                continue
            if rt_range is not None:
                if spec.rt < rt_range[0] or spec.rt > rt_range[1]:
                    parse_scans(scan)
                    continue

            exp.add_spectrum(spec)
            loaded += 1

            if progress_callback:
                if not progress_callback(loaded, -1):
                    return

            # Parse nested scans (MS/MS)
            parse_scans(scan)

    parse_scans(ms_run)
    return exp


def _parse_mzxml_scan(elem: ET.Element, ns: Dict[str, str]) -> Spectrum:
    """Parse an mzXML scan element."""
    spec = Spectrum()

    spec.ms_level = int(elem.get("msLevel", 1))

    # Parse RT
    rt_str = elem.get("retentionTime", "0")
    if rt_str.startswith("PT") and rt_str.endswith("S"):
        spec.rt = float(rt_str[2:-1])
    elif rt_str.startswith("PT") and rt_str.endswith("M"):
        spec.rt = float(rt_str[2:-1]) * 60.0
    else:
        spec.rt = float(rt_str)

    # Polarity
    polarity = elem.get("polarity", "")
    if polarity in ("+", "positive"):
        spec.polarity = Polarity.POSITIVE
    elif polarity in ("-", "negative"):
        spec.polarity = Polarity.NEGATIVE

    # Centroided
    if elem.get("centroided", "0") in ("1", "true"):
        spec.spectrum_type = SpectrumType.CENTROID
    else:
        spec.spectrum_type = SpectrumType.PROFILE

    # Parse precursors
    for prec_elem in elem.findall("mzxml:precursorMz", ns) or elem.findall("precursorMz"):
        prec = Precursor()
        prec.mz = float(prec_elem.text or 0)
        prec.intensity = float(prec_elem.get("precursorIntensity", 0))
        prec.charge = int(prec_elem.get("precursorCharge", 0))
        spec.precursors.append(prec)

    # Parse peaks
    peaks_elem = elem.find("mzxml:peaks", ns) or elem.find("peaks")
    if peaks_elem is not None and peaks_elem.text:
        precision = int(peaks_elem.get("precision", 32))
        byte_order = peaks_elem.get("byteOrder", "network")
        compression = peaks_elem.get("compressionType", "")
        pair_order = peaks_elem.get("pairOrder", "m/z-int")

        is_64bit = (precision == 64)
        is_little_endian = (byte_order != "network")
        is_compressed = (compression == "zlib")

        values = decode_binary(peaks_elem.text, is_64bit, is_compressed, is_little_endian)

        # Values are interleaved
        mz_data = []
        intensity_data = []
        mz_first = (pair_order != "int-mz")

        for i in range(0, len(values) - 1, 2):
            if mz_first:
                mz_data.append(values[i])
                intensity_data.append(values[i + 1])
            else:
                intensity_data.append(values[i])
                mz_data.append(values[i + 1])

        spec.mz = mz_data
        spec.intensity = intensity_data

    return spec


def save_mztab(
    features: Any,
    filename: str,
    study_id: str = "lcms_study",
    run_id: str = "run_1",
) -> None:
    """
    Save features to mzTab format.

    Args:
        features: FeatureMap or similar with features
        filename: Output filename
        study_id: Study identifier
        run_id: MS run identifier
    """
    from .feature import FeatureMap

    with open(filename, "w") as f:
        # Metadata section
        f.write("MTD\tmzTab-version\t1.0.0\n")
        f.write("MTD\tmzTab-mode\tSummary\n")
        f.write("MTD\tmzTab-type\tQuantification\n")
        f.write(f"MTD\tdescription\t{study_id}\n")
        f.write(f"MTD\tms_run[1]-location\t{features.source_file if hasattr(features, 'source_file') else 'unknown'}\n")

        # Small molecule header
        f.write("\n")
        f.write("SMH\tSML_ID\tidentifier\tsmiles\tmass\tcharge\tretention_time\topt_peak_intensity\topt_peak_area\n")

        # Small molecule data
        for i, feat in enumerate(features):
            f.write(f"SML\t{i+1}\t")
            f.write(f"feature_{i+1}\t")
            f.write("\t")  # SMILES (empty)
            f.write(f"{feat.mz:.6f}\t")
            f.write(f"{feat.charge if feat.charge else ''}\t")
            f.write(f"{feat.rt:.2f}\t")
            f.write(f"{feat.intensity:.0f}\t")
            f.write(f"{feat.volume:.0f}\n")

    print(f"Saved {len(list(features))} features to {filename}")
