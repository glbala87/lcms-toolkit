package org.lcms.core;

import java.util.*;

/**
 * Matrix of features (rows) x samples (columns) with intensities.
 * Represents aligned features across multiple LC-MS runs.
 */
public class ConsensusMap {

    private final List<String> featureIds;
    private final List<String> sampleNames;
    private double[][] intensityMatrix;

    public ConsensusMap(List<String> featureIds, List<String> sampleNames,
                        double[][] intensityMatrix) {
        this.featureIds = new ArrayList<>(featureIds);
        this.sampleNames = new ArrayList<>(sampleNames);
        this.intensityMatrix = intensityMatrix;
    }

    public int getNumFeatures() { return featureIds.size(); }
    public int getNumSamples() { return sampleNames.size(); }
    public List<String> getFeatureIds() { return featureIds; }
    public List<String> getSampleNames() { return sampleNames; }

    public double getIntensity(int feature, int sample) {
        return intensityMatrix[feature][sample];
    }

    public void setIntensity(int feature, int sample, double value) {
        intensityMatrix[feature][sample] = value;
    }

    public double[][] getIntensityMatrix() { return intensityMatrix; }

    /**
     * Filter features by minimum presence across samples.
     */
    public ConsensusMap filterByPresence(double minFraction) {
        List<String> filteredIds = new ArrayList<>();
        List<double[]> filteredRows = new ArrayList<>();

        for (int i = 0; i < getNumFeatures(); i++) {
            int present = 0;
            for (int j = 0; j < getNumSamples(); j++) {
                if (intensityMatrix[i][j] > 0) present++;
            }
            double fraction = (double) present / getNumSamples();
            if (fraction >= minFraction) {
                filteredIds.add(featureIds.get(i));
                filteredRows.add(intensityMatrix[i].clone());
            }
        }

        double[][] newMatrix = filteredRows.toArray(new double[0][]);
        return new ConsensusMap(filteredIds, sampleNames, newMatrix);
    }

    /**
     * Fill missing (zero) values with column minimum.
     */
    public ConsensusMap fillMissing() {
        double[][] filled = new double[getNumFeatures()][getNumSamples()];

        for (int j = 0; j < getNumSamples(); j++) {
            double colMin = Double.MAX_VALUE;
            for (int i = 0; i < getNumFeatures(); i++) {
                if (intensityMatrix[i][j] > 0 && intensityMatrix[i][j] < colMin) {
                    colMin = intensityMatrix[i][j];
                }
            }
            if (colMin == Double.MAX_VALUE) colMin = 1.0;

            for (int i = 0; i < getNumFeatures(); i++) {
                filled[i][j] = intensityMatrix[i][j] > 0
                        ? intensityMatrix[i][j] : colMin;
            }
        }

        return new ConsensusMap(featureIds, sampleNames, filled);
    }

    /**
     * Apply median normalization.
     */
    public ConsensusMap medianNormalize() {
        double[][] normalized = Quantification.medianNormalization(intensityMatrix);
        return new ConsensusMap(featureIds, sampleNames, normalized);
    }

    /**
     * Log2-transform intensities.
     */
    public ConsensusMap log2Transform() {
        double[][] transformed = new double[getNumFeatures()][getNumSamples()];
        for (int i = 0; i < getNumFeatures(); i++) {
            for (int j = 0; j < getNumSamples(); j++) {
                transformed[i][j] = Math.log(intensityMatrix[i][j] + 1) / Math.log(2);
            }
        }
        return new ConsensusMap(featureIds, sampleNames, transformed);
    }

    @Override
    public String toString() {
        return String.format("ConsensusMap(%d features x %d samples)",
                getNumFeatures(), getNumSamples());
    }
}
