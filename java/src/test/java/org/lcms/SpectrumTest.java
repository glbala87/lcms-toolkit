package org.lcms;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.BeforeEach;
import static org.junit.jupiter.api.Assertions.*;

import org.lcms.core.Spectrum;

/**
 * JUnit 5 tests for core Spectrum class.
 */
class SpectrumTest {

    @Test
    void testCreateEmpty() {
        Spectrum spec = new Spectrum();
        assertEquals(0, spec.getMzArray().length);
        assertEquals(0, spec.getIntensityArray().length);
    }

    @Test
    void testCreateWithData() {
        double[] mz = {100.0, 200.0, 300.0};
        double[] ints = {1000.0, 2000.0, 500.0};
        Spectrum spec = new Spectrum(mz, ints);
        assertEquals(3, spec.getMzArray().length);
        assertArrayEquals(mz, spec.getMzArray(), 1e-10);
    }

    @Test
    void testBasePeakMz() {
        double[] mz = {100.0, 200.0, 300.0};
        double[] ints = {1000.0, 2000.0, 500.0};
        Spectrum spec = new Spectrum(mz, ints);
        assertEquals(200.0, spec.getBasePeakMz(), 1e-10);
    }

    @Test
    void testTotalIonCurrent() {
        double[] mz = {100.0, 200.0, 300.0};
        double[] ints = {1000.0, 2000.0, 500.0};
        Spectrum spec = new Spectrum(mz, ints);
        assertEquals(3500.0, spec.getTotalIonCurrent(), 1e-10);
    }
}
