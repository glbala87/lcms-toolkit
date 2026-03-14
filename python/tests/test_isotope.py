"""Tests for isotope detection and deconvolution."""

import pytest
import numpy as np

from pylcms.isotope import (
    IsotopePattern, DeconvolutedMass,
    detect_isotope_patterns, deconvolute_spectrum,
    averagine_distribution, assign_charge_state,
)
from pylcms.spectrum import Spectrum


class TestAveragine:
    def test_averagine_distribution(self):
        dist = averagine_distribution(1000.0, n_peaks=5)
        assert len(dist) == 5
        assert abs(sum(dist) - 1.0) < 0.01
        assert dist[0] > dist[-1]  # Monoisotopic should be most intense for small masses

    def test_averagine_large_mass(self):
        dist = averagine_distribution(5000.0, n_peaks=8)
        assert len(dist) == 8
        assert sum(dist) > 0.95


class TestIsotopeDetection:
    def _make_isotope_spectrum(self, mono_mz=500.0, charge=2, n_isotopes=5):
        """Create a synthetic spectrum with an isotope pattern."""
        spacing = 1.003355 / charge
        mz_list = []
        int_list = []
        dist = averagine_distribution(mono_mz * charge, n_peaks=n_isotopes)
        for i in range(n_isotopes):
            mz_list.append(mono_mz + i * spacing)
            int_list.append(dist[i] * 10000)
        # Add noise peaks
        for _ in range(50):
            mz_list.append(np.random.uniform(400, 600))
            int_list.append(np.random.uniform(10, 100))

        idx = np.argsort(mz_list)
        return Spectrum(
            mz_array=np.array(mz_list)[idx],
            intensity_array=np.array(int_list)[idx],
        )

    def test_detect_patterns(self):
        spec = self._make_isotope_spectrum(mono_mz=500.0, charge=2)
        patterns = detect_isotope_patterns(spec)
        assert len(patterns) >= 1

    def test_assign_charge(self):
        charge = assign_charge_state(
            [500.0, 500.5016, 501.0032],
            [10000, 8000, 4000],
        )
        assert charge == 2


class TestDeconvolution:
    def test_deconvolute(self):
        mono_mz = 500.0
        charge = 2
        spacing = 1.003355 / charge
        mz = np.array([mono_mz + i * spacing for i in range(5)])
        ints = np.array([10000, 8000, 4000, 1500, 500])

        # Add noise
        noise_mz = np.random.uniform(400, 600, 50)
        noise_int = np.random.uniform(10, 100, 50)
        all_mz = np.concatenate([mz, noise_mz])
        all_int = np.concatenate([ints, noise_int])
        idx = np.argsort(all_mz)

        spec = Spectrum(mz_array=all_mz[idx], intensity_array=all_int[idx])
        results = deconvolute_spectrum(spec)
        assert len(results) >= 1
