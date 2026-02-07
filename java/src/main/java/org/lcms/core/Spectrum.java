package org.lcms.core;

import java.util.*;

/**
 * Represents a mass spectrum with m/z and intensity data.
 */
public class Spectrum {
    private double[] mz;
    private double[] intensity;
    private int index;
    private String nativeId;
    private int msLevel;
    private double retentionTime;
    private SpectrumType type;
    private Polarity polarity;
    private List<Precursor> precursors;
    private Map<String, String> metadata;

    // Cached statistics
    private double tic;
    private double basePeakIntensity;
    private double basePeakMz;
    private double mzMin;
    private double mzMax;

    public Spectrum() {
        this.mz = new double[0];
        this.intensity = new double[0];
        this.msLevel = 1;
        this.retentionTime = 0.0;
        this.type = SpectrumType.UNKNOWN;
        this.polarity = Polarity.UNKNOWN;
        this.precursors = new ArrayList<>();
        this.metadata = new HashMap<>();
    }

    public Spectrum(double[] mz, double[] intensity) {
        this();
        setData(mz, intensity);
    }

    public void setData(double[] mz, double[] intensity) {
        if (mz.length != intensity.length) {
            throw new IllegalArgumentException("m/z and intensity arrays must have same length");
        }
        this.mz = mz.clone();
        this.intensity = intensity.clone();
        updateStatistics();
    }

    private void updateStatistics() {
        if (mz.length == 0) {
            tic = 0;
            basePeakIntensity = 0;
            basePeakMz = 0;
            mzMin = 0;
            mzMax = 0;
            return;
        }

        tic = 0;
        basePeakIntensity = Double.NEGATIVE_INFINITY;
        mzMin = Double.POSITIVE_INFINITY;
        mzMax = Double.NEGATIVE_INFINITY;

        for (int i = 0; i < mz.length; i++) {
            tic += intensity[i];
            if (intensity[i] > basePeakIntensity) {
                basePeakIntensity = intensity[i];
                basePeakMz = mz[i];
            }
            if (mz[i] < mzMin) mzMin = mz[i];
            if (mz[i] > mzMax) mzMax = mz[i];
        }
    }

    // Getters and setters
    public int size() { return mz.length; }
    public boolean isEmpty() { return mz.length == 0; }

    public double[] getMz() { return mz; }
    public double[] getIntensity() { return intensity; }

    public double getMzAt(int i) { return mz[i]; }
    public double getIntensityAt(int i) { return intensity[i]; }

    public int getIndex() { return index; }
    public void setIndex(int index) { this.index = index; }

    public String getNativeId() { return nativeId; }
    public void setNativeId(String nativeId) { this.nativeId = nativeId; }

    public int getMsLevel() { return msLevel; }
    public void setMsLevel(int msLevel) { this.msLevel = msLevel; }

    public double getRetentionTime() { return retentionTime; }
    public void setRetentionTime(double retentionTime) { this.retentionTime = retentionTime; }

    public SpectrumType getType() { return type; }
    public void setType(SpectrumType type) { this.type = type; }

    public Polarity getPolarity() { return polarity; }
    public void setPolarity(Polarity polarity) { this.polarity = polarity; }

    public List<Precursor> getPrecursors() { return precursors; }
    public void addPrecursor(Precursor p) { precursors.add(p); }

    public Map<String, String> getMetadata() { return metadata; }

    public double getTic() { return tic; }
    public double getBasePeakIntensity() { return basePeakIntensity; }
    public double getBasePeakMz() { return basePeakMz; }
    public double getMzMin() { return mzMin; }
    public double getMzMax() { return mzMax; }

    /**
     * Check if spectrum is sorted by m/z.
     */
    public boolean isSorted() {
        for (int i = 1; i < mz.length; i++) {
            if (mz[i] < mz[i-1]) return false;
        }
        return true;
    }

    /**
     * Sort spectrum by m/z.
     */
    public void sortByMz() {
        if (mz.length <= 1) return;

        Integer[] indices = new Integer[mz.length];
        for (int i = 0; i < indices.length; i++) indices[i] = i;

        Arrays.sort(indices, Comparator.comparingDouble(i -> mz[i]));

        double[] newMz = new double[mz.length];
        double[] newIntensity = new double[intensity.length];
        for (int i = 0; i < indices.length; i++) {
            newMz[i] = mz[indices[i]];
            newIntensity[i] = intensity[indices[i]];
        }
        mz = newMz;
        intensity = newIntensity;
    }

    /**
     * Find index of peak closest to target m/z.
     */
    public int findNearestMz(double target) {
        if (mz.length == 0) {
            throw new IllegalStateException("Cannot find peak in empty spectrum");
        }

        int nearest = 0;
        double minDist = Math.abs(mz[0] - target);

        for (int i = 1; i < mz.length; i++) {
            double dist = Math.abs(mz[i] - target);
            if (dist < minDist) {
                minDist = dist;
                nearest = i;
            }
        }
        return nearest;
    }

    /**
     * Extract peaks within m/z range.
     */
    public Spectrum extractRange(double mzLow, double mzHigh) {
        List<Double> newMz = new ArrayList<>();
        List<Double> newIntensity = new ArrayList<>();

        for (int i = 0; i < mz.length; i++) {
            if (mz[i] >= mzLow && mz[i] <= mzHigh) {
                newMz.add(mz[i]);
                newIntensity.add(intensity[i]);
            }
        }

        Spectrum result = new Spectrum();
        result.setMsLevel(msLevel);
        result.setRetentionTime(retentionTime);
        result.setType(type);
        result.setPolarity(polarity);
        result.setData(
            newMz.stream().mapToDouble(Double::doubleValue).toArray(),
            newIntensity.stream().mapToDouble(Double::doubleValue).toArray()
        );
        return result;
    }

    @Override
    public String toString() {
        return String.format("Spectrum(size=%d, msLevel=%d, rt=%.2fs, mz=[%.2f, %.2f])",
            size(), msLevel, retentionTime, mzMin, mzMax);
    }
}
