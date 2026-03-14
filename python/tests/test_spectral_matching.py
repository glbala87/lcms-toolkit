"""Tests for spectral matching and similarity scoring."""

import pytest
import numpy as np

from pylcms.spectral_matching import (
    cosine_similarity,
    modified_cosine_similarity,
    spectral_entropy_similarity,
    SpectralLibrary,
    SpectralMatch,
)
from pylcms.spectrum import Spectrum


def _make_spectrum(mz_values, intensity_values):
    return Spectrum(
        mz_array=np.array(mz_values, dtype=float),
        intensity_array=np.array(intensity_values, dtype=float),
    )


class TestCosineSimilarity:
    def test_identical_spectra(self):
        spec = _make_spectrum([100, 200, 300], [1000, 2000, 500])
        score = cosine_similarity(spec, spec)
        assert abs(score - 1.0) < 0.01

    def test_orthogonal_spectra(self):
        spec1 = _make_spectrum([100, 200], [1000, 0])
        spec2 = _make_spectrum([300, 400], [0, 1000])
        score = cosine_similarity(spec1, spec2)
        assert score < 0.1

    def test_similar_spectra(self):
        spec1 = _make_spectrum([100, 200, 300], [1000, 2000, 500])
        spec2 = _make_spectrum([100.001, 200.001, 300.001], [900, 2100, 480])
        score = cosine_similarity(spec1, spec2, tolerance=0.01)
        assert score > 0.9


class TestModifiedCosine:
    def test_with_precursor_shift(self):
        spec1 = _make_spectrum([100, 200, 300], [1000, 2000, 500])
        spec2 = _make_spectrum([114, 214, 314], [1000, 2000, 500])
        score = modified_cosine_similarity(
            spec1, spec2,
            precursor_mz1=400.0, precursor_mz2=414.0,
            tolerance=0.02,
        )
        assert score > 0.8


class TestEntropySimiliarity:
    def test_identical(self):
        spec = _make_spectrum([100, 200, 300], [1000, 2000, 500])
        score = spectral_entropy_similarity(spec, spec)
        assert score > 0.95

    def test_different(self):
        spec1 = _make_spectrum([100, 200], [1000, 500])
        spec2 = _make_spectrum([300, 400], [800, 1200])
        score = spectral_entropy_similarity(spec1, spec2)
        assert score < 0.2


class TestSpectralLibrary:
    def test_create_library(self):
        lib = SpectralLibrary()
        assert len(lib) == 0

    def test_add_and_search(self):
        lib = SpectralLibrary()
        ref = _make_spectrum([100, 200, 300], [1000, 2000, 500])
        lib.add_spectrum(ref, name="test_compound", precursor_mz=400.0)
        assert len(lib) == 1

        query = _make_spectrum([100.001, 200.001, 300.001], [950, 2050, 520])
        matches = lib.search(query, precursor_mz=400.0, top_k=5)
        assert len(matches) >= 1
        assert matches[0].score > 0.8
