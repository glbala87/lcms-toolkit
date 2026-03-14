"""Tests for Spectrum and related data structures."""

import pytest
import numpy as np

from pylcms.spectrum import Spectrum, SpectrumType, Polarity
from pylcms.peak import Peak, PeakList


class TestSpectrum:
    def test_create_empty(self):
        spec = Spectrum()
        assert len(spec.mz_array) == 0
        assert len(spec.intensity_array) == 0

    def test_create_with_data(self):
        mz = np.array([100.0, 200.0, 300.0])
        ints = np.array([1000.0, 2000.0, 500.0])
        spec = Spectrum(mz_array=mz, intensity_array=ints)
        assert len(spec.mz_array) == 3
        np.testing.assert_array_equal(spec.mz_array, mz)

    def test_spectrum_type(self):
        assert SpectrumType.MS1.value in ("ms1", "MS1", 1)
        assert SpectrumType.MS2.value in ("ms2", "MS2", 2)

    def test_polarity(self):
        assert Polarity.POSITIVE is not None
        assert Polarity.NEGATIVE is not None


class TestPeak:
    def test_create_peak(self):
        peak = Peak(mz=500.5, intensity=1234.0, rt=120.0)
        assert peak.mz == 500.5
        assert peak.intensity == 1234.0
        assert peak.rt == 120.0

    def test_peak_list(self):
        peaks = PeakList()
        peaks.add(Peak(mz=100.0, intensity=500.0))
        peaks.add(Peak(mz=200.0, intensity=1000.0))
        assert len(peaks) == 2


class TestAlgorithms:
    def test_pick_peaks(self):
        from pylcms.algorithms import pick_peaks

        mz = np.linspace(100, 200, 1000)
        ints = np.zeros(1000)
        # Add a Gaussian peak
        center = 500
        for i in range(1000):
            ints[i] = 1000 * np.exp(-0.5 * ((i - center) / 10) ** 2)

        spec = Spectrum(mz_array=mz, intensity_array=ints)
        peaks = pick_peaks(spec)
        assert len(peaks) >= 1

    def test_smooth_spectrum(self):
        from pylcms.algorithms import smooth_spectrum

        mz = np.linspace(100, 200, 100)
        ints = np.random.rand(100) * 1000
        spec = Spectrum(mz_array=mz, intensity_array=ints)
        smoothed = smooth_spectrum(spec)
        assert len(smoothed.mz_array) == len(mz)

    def test_estimate_noise(self):
        from pylcms.algorithms import estimate_noise

        mz = np.linspace(100, 200, 1000)
        noise = np.random.normal(0, 10, 1000)
        spec = Spectrum(mz_array=mz, intensity_array=np.abs(noise))
        noise_level = estimate_noise(spec)
        assert noise_level > 0
