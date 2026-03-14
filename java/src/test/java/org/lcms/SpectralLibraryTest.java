package org.lcms;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.BeforeEach;
import static org.junit.jupiter.api.Assertions.*;

import org.lcms.core.SpectralLibrary;
import org.lcms.core.SpectralMatch;
import org.lcms.core.Spectrum;

import java.util.List;

/**
 * JUnit 5 tests for SpectralLibrary.
 */
class SpectralLibraryTest {

    private SpectralLibrary library;

    @BeforeEach
    void setUp() {
        library = new SpectralLibrary();
    }

    @Test
    void testEmptyLibrary() {
        assertEquals(0, library.size());
    }

    @Test
    void testAddSpectrum() {
        double[] mz = {100.0, 200.0, 300.0};
        double[] ints = {1000.0, 2000.0, 500.0};
        Spectrum spec = new Spectrum(mz, ints);
        library.addSpectrum("compound1", 400.0, spec);
        assertEquals(1, library.size());
    }

    @Test
    void testSearchIdentical() {
        double[] mz = {100.0, 200.0, 300.0};
        double[] ints = {1000.0, 2000.0, 500.0};
        Spectrum ref = new Spectrum(mz, ints);
        library.addSpectrum("compound1", 400.0, ref);

        Spectrum query = new Spectrum(mz, ints);
        List<SpectralMatch> matches = library.search(query, 400.0, 5);
        assertFalse(matches.isEmpty());
        assertTrue(matches.get(0).getScore() > 0.9);
    }

    @Test
    void testSearchNoMatch() {
        double[] refMz = {100.0, 200.0};
        double[] refInts = {1000.0, 2000.0};
        library.addSpectrum("compound1", 400.0, new Spectrum(refMz, refInts));

        double[] queryMz = {500.0, 600.0};
        double[] queryInts = {1000.0, 2000.0};
        List<SpectralMatch> matches = library.search(
            new Spectrum(queryMz, queryInts), 700.0, 5);
        // Should return empty or low-scoring matches
        assertTrue(matches.isEmpty() || matches.get(0).getScore() < 0.3);
    }
}
