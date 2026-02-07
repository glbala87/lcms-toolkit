package org.lcms.core;

import java.util.*;

/**
 * Represents a chromatogram (intensity vs retention time).
 */
public class Chromatogram {
    private double[] rt;
    private double[] intensity;
    private int index;
    private String nativeId;
    private ChromatogramType type;
    private Polarity polarity;
    private double targetMz;
    private double mzTolerance;
    private double precursorMz;
    private double productMz;
    private Map<String, String> metadata;

    // Cached statistics
    private double maxIntensity;
    private double apexRt;
    private double rtMin;
    private double rtMax;

    public Chromatogram() {
        this.rt = new double[0];
        this.intensity = new double[0];
        this.type = ChromatogramType.UNKNOWN;
        this.polarity = Polarity.UNKNOWN;
        this.metadata = new HashMap<>();
    }

    public Chromatogram(double[] rt, double[] intensity) {
        this();
        setData(rt, intensity);
    }

    public void setData(double[] rt, double[] intensity) {
        if (rt.length != intensity.length) {
            throw new IllegalArgumentException("RT and intensity arrays must have same length");
        }
        this.rt = rt.clone();
        this.intensity = intensity.clone();
        updateStatistics();
    }

    private void updateStatistics() {
        if (rt.length == 0) {
            maxIntensity = 0;
            apexRt = 0;
            rtMin = 0;
            rtMax = 0;
            return;
        }

        maxIntensity = Double.NEGATIVE_INFINITY;
        rtMin = Double.POSITIVE_INFINITY;
        rtMax = Double.NEGATIVE_INFINITY;

        for (int i = 0; i < rt.length; i++) {
            if (intensity[i] > maxIntensity) {
                maxIntensity = intensity[i];
                apexRt = rt[i];
            }
            if (rt[i] < rtMin) rtMin = rt[i];
            if (rt[i] > rtMax) rtMax = rt[i];
        }
    }

    // Getters and setters
    public int size() { return rt.length; }
    public boolean isEmpty() { return rt.length == 0; }

    public double[] getRt() { return rt; }
    public double[] getIntensity() { return intensity; }

    public double getRtAt(int i) { return rt[i]; }
    public double getIntensityAt(int i) { return intensity[i]; }

    public int getIndex() { return index; }
    public void setIndex(int index) { this.index = index; }

    public String getNativeId() { return nativeId; }
    public void setNativeId(String nativeId) { this.nativeId = nativeId; }

    public ChromatogramType getType() { return type; }
    public void setType(ChromatogramType type) { this.type = type; }

    public Polarity getPolarity() { return polarity; }
    public void setPolarity(Polarity polarity) { this.polarity = polarity; }

    public double getTargetMz() { return targetMz; }
    public void setTargetMz(double targetMz) { this.targetMz = targetMz; }

    public double getMzTolerance() { return mzTolerance; }
    public void setMzTolerance(double tolerance) { this.mzTolerance = tolerance; }

    public double getPrecursorMz() { return precursorMz; }
    public void setPrecursorMz(double mz) { this.precursorMz = mz; }

    public double getProductMz() { return productMz; }
    public void setProductMz(double mz) { this.productMz = mz; }

    public Map<String, String> getMetadata() { return metadata; }

    public double getMaxIntensity() { return maxIntensity; }
    public double getApexRt() { return apexRt; }
    public double getRtMin() { return rtMin; }
    public double getRtMax() { return rtMax; }

    /**
     * Check if chromatogram is sorted by RT.
     */
    public boolean isSorted() {
        for (int i = 1; i < rt.length; i++) {
            if (rt[i] < rt[i-1]) return false;
        }
        return true;
    }

    /**
     * Sort chromatogram by RT.
     */
    public void sortByRt() {
        if (rt.length <= 1) return;

        Integer[] indices = new Integer[rt.length];
        for (int i = 0; i < indices.length; i++) indices[i] = i;

        Arrays.sort(indices, Comparator.comparingDouble(i -> rt[i]));

        double[] newRt = new double[rt.length];
        double[] newIntensity = new double[intensity.length];
        for (int i = 0; i < indices.length; i++) {
            newRt[i] = rt[indices[i]];
            newIntensity[i] = intensity[indices[i]];
        }
        rt = newRt;
        intensity = newIntensity;
    }

    /**
     * Compute area under the curve using trapezoidal integration.
     */
    public double computeArea() {
        if (rt.length < 2) return 0;

        double area = 0;
        for (int i = 1; i < rt.length; i++) {
            double dt = rt[i] - rt[i-1];
            area += 0.5 * (intensity[i] + intensity[i-1]) * dt;
        }
        return area;
    }

    /**
     * Extract points within RT range.
     */
    public Chromatogram extractRange(double rtLow, double rtHigh) {
        List<Double> newRt = new ArrayList<>();
        List<Double> newIntensity = new ArrayList<>();

        for (int i = 0; i < rt.length; i++) {
            if (rt[i] >= rtLow && rt[i] <= rtHigh) {
                newRt.add(rt[i]);
                newIntensity.add(intensity[i]);
            }
        }

        Chromatogram result = new Chromatogram();
        result.setType(type);
        result.setPolarity(polarity);
        result.setTargetMz(targetMz);
        result.setData(
            newRt.stream().mapToDouble(Double::doubleValue).toArray(),
            newIntensity.stream().mapToDouble(Double::doubleValue).toArray()
        );
        return result;
    }

    @Override
    public String toString() {
        return String.format("Chromatogram(size=%d, type=%s, rt=[%.1f, %.1f])",
            size(), type, rtMin, rtMax);
    }
}
