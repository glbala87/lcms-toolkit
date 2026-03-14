"""Tests for peptide identification module."""

import pytest
import numpy as np

from pylcms.identification import (
    calculate_peptide_mass,
    generate_theoretical_fragments,
    target_decoy_fdr,
    PeptideSpectrumMatch,
)


class TestPeptideMass:
    def test_known_mass(self):
        # Glycine (G) monoisotopic mass = 57.021464
        # Water = 18.010565
        mass = calculate_peptide_mass("G")
        assert abs(mass - (57.021464 + 18.010565)) < 0.01

    def test_longer_peptide(self):
        mass = calculate_peptide_mass("PEPTIDE")
        assert mass > 700  # Should be around 799.36

    def test_empty(self):
        mass = calculate_peptide_mass("")
        assert abs(mass - 18.010565) < 0.01  # Just water


class TestFragments:
    def test_generate_fragments(self):
        fragments = generate_theoretical_fragments("PEPTIDE", charge=1)
        assert len(fragments) > 0

        # Should have b and y ions
        ion_types = set(f[1] for f in fragments)
        assert "b" in ion_types or any("b" in t for t in ion_types)

    def test_fragment_count(self):
        # For a 7-residue peptide, expect 6 b-ions and 6 y-ions at minimum
        fragments = generate_theoretical_fragments("PEPTIDE", charge=1)
        assert len(fragments) >= 12


class TestFDR:
    def test_target_decoy(self):
        # Create mock PSMs
        psms = []
        # 10 good target hits
        for i in range(10):
            psms.append(PeptideSpectrumMatch(
                sequence=f"PEPTIDE{i}",
                score=50.0 + i,
                is_decoy=False,
            ))
        # 2 decoy hits
        for i in range(2):
            psms.append(PeptideSpectrumMatch(
                sequence=f"DECOY{i}",
                score=20.0 + i,
                is_decoy=True,
            ))

        filtered = target_decoy_fdr(psms, fdr_threshold=0.05)
        assert len(filtered) <= len(psms)
        # All filtered should be targets at this FDR
        assert all(not p.is_decoy for p in filtered)
