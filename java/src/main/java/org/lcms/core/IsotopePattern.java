package org.lcms.core;

import java.util.*;

/**
 * Represents a detected isotope envelope in a mass spectrum.
 */
public class IsotopePattern {

    public static final double NEUTRON_MASS = 1.003355;
    public static final double PROTON_MASS = 1.007276;
    public static final double AVERAGINE_MASS = 111.1254;

    private double monoisotopicMz;
    private int charge;
    private List<double[]> peaks; // Each element is {mz, intensity}
    private double score;
    private double retentionTime;

    public IsotopePattern() {
        this.charge = 1;
        this.peaks = new ArrayList<>();
    }

    public IsotopePattern(double monoisotopicMz, int charge, List<double[]> peaks, double score) {
        this.monoisotopicMz = monoisotopicMz;
        this.charge = charge;
        this.peaks = peaks != null ? new ArrayList<>(peaks) : new ArrayList<>();
        this.score = score;
    }

    public double getMonoisotopicMz() { return monoisotopicMz; }
    public void setMonoisotopicMz(double mz) { this.monoisotopicMz = mz; }

    public int getCharge() { return charge; }
    public void setCharge(int charge) { this.charge = charge; }

    public List<double[]> getPeaks() { return peaks; }
    public void setPeaks(List<double[]> peaks) { this.peaks = peaks; }

    public double getScore() { return score; }
    public void setScore(double score) { this.score = score; }

    public double getRetentionTime() { return retentionTime; }
    public void setRetentionTime(double rt) { this.retentionTime = rt; }

    public int getNumPeaks() { return peaks.size(); }

    /**
     * Calculate neutral mass from monoisotopic m/z and charge.
     */
    public double getNeutralMass() {
        return (monoisotopicMz - PROTON_MASS) * Math.abs(charge);
    }

    /**
     * Total intensity across all isotope peaks.
     */
    public double getTotalIntensity() {
        double total = 0.0;
        for (double[] peak : peaks) {
            total += peak[1];
        }
        return total;
    }

    /**
     * Generate theoretical isotope distribution using averagine model.
     */
    public static double[] averagineDistribution(double mass, int numPeaks) {
        double[] distribution = new double[numPeaks];
        if (numPeaks <= 0) return distribution;

        if (mass <= 0) {
            distribution[0] = 1.0;
            return distribution;
        }

        double nUnits = mass / AVERAGINE_MASS;
        double lambda = nUnits * (
            4.9384 * 0.0111 +
            1.3577 * 0.00366 +
            1.4773 * 0.00205 +
            0.0417 * 0.0425 +
            7.7583 * 0.00012
        );

        for (int i = 0; i < numPeaks; i++) {
            double logProb = -lambda + i * Math.log(lambda);
            double logFact = 0.0;
            for (int j = 2; j <= i; j++) {
                logFact += Math.log(j);
            }
            distribution[i] = Math.exp(logProb - logFact);
        }

        // Normalize to max = 1
        double max = 0.0;
        for (double v : distribution) {
            if (v > max) max = v;
        }
        if (max > 0) {
            for (int i = 0; i < distribution.length; i++) {
                distribution[i] /= max;
            }
        }

        return distribution;
    }

    /**
     * Compute cosine similarity between two vectors.
     */
    public static double cosineSimilarity(double[] a, double[] b) {
        int n = Math.min(a.length, b.length);
        if (n == 0) return 0.0;

        double dot = 0, normA = 0, normB = 0;
        for (int i = 0; i < n; i++) {
            dot += a[i] * b[i];
            normA += a[i] * a[i];
            normB += b[i] * b[i];
        }

        normA = Math.sqrt(normA);
        normB = Math.sqrt(normB);
        if (normA == 0 || normB == 0) return 0.0;
        return dot / (normA * normB);
    }

    @Override
    public String toString() {
        return String.format("IsotopePattern(mz=%.4f, z=%d, peaks=%d, score=%.3f)",
                monoisotopicMz, charge, peaks.size(), score);
    }
}
