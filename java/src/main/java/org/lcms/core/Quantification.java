package org.lcms.core;

import java.util.*;

/**
 * Label-free quantification: feature alignment, normalization, and differential analysis.
 */
public class Quantification {

    private double mzTolerance = 0.01;
    private double rtTolerance = 30.0;

    public Quantification() {}

    public Quantification(double mzTolerance, double rtTolerance) {
        this.mzTolerance = mzTolerance;
        this.rtTolerance = rtTolerance;
    }

    public void setMzTolerance(double v) { this.mzTolerance = v; }
    public void setRtTolerance(double v) { this.rtTolerance = v; }

    /**
     * Align features across multiple runs into a consensus map.
     */
    public ConsensusMap alignFeatures(List<double[][]> featureMaps,
                                       List<String> sampleNames) {
        // featureMaps: each is double[n][3] with {mz, rt, intensity}
        int nSamples = featureMaps.size();
        if (nSamples == 0) {
            return new ConsensusMap(new ArrayList<>(), sampleNames, new double[0][0]);
        }

        List<double[]> consensusMz = new ArrayList<>();
        List<double[]> consensusRt = new ArrayList<>();
        List<Map<Integer, Double>> consensusIntensities = new ArrayList<>();

        // Seed with first sample
        double[][] first = featureMaps.get(0);
        for (double[] feat : first) {
            consensusMz.add(new double[]{feat[0]});
            consensusRt.add(new double[]{feat[1]});
            Map<Integer, Double> intensities = new HashMap<>();
            intensities.put(0, feat[2]);
            consensusIntensities.add(intensities);
        }

        // Match features from remaining samples
        for (int s = 1; s < nSamples; s++) {
            double[][] features = featureMaps.get(s);
            Set<Integer> usedConsensus = new HashSet<>();

            for (double[] feat : features) {
                int bestMatch = -1;
                double bestDist = Double.MAX_VALUE;

                for (int ci = 0; ci < consensusMz.size(); ci++) {
                    if (usedConsensus.contains(ci)) continue;

                    double mzDiff = Math.abs(feat[0] - consensusMz.get(ci)[0]);
                    double rtDiff = Math.abs(feat[1] - consensusRt.get(ci)[0]);

                    if (mzDiff <= mzTolerance && rtDiff <= rtTolerance) {
                        double dist = mzDiff / mzTolerance + rtDiff / rtTolerance;
                        if (dist < bestDist) {
                            bestDist = dist;
                            bestMatch = ci;
                        }
                    }
                }

                if (bestMatch >= 0) {
                    consensusIntensities.get(bestMatch).put(s, feat[2]);
                    usedConsensus.add(bestMatch);
                } else {
                    consensusMz.add(new double[]{feat[0]});
                    consensusRt.add(new double[]{feat[1]});
                    Map<Integer, Double> intensities = new HashMap<>();
                    intensities.put(s, feat[2]);
                    consensusIntensities.add(intensities);
                }
            }
        }

        // Build matrix
        int nFeatures = consensusMz.size();
        double[][] matrix = new double[nFeatures][nSamples];
        List<String> featureIds = new ArrayList<>();

        for (int i = 0; i < nFeatures; i++) {
            featureIds.add(String.format("F%06d", i));
            Map<Integer, Double> intensities = consensusIntensities.get(i);
            for (Map.Entry<Integer, Double> entry : intensities.entrySet()) {
                matrix[i][entry.getKey()] = entry.getValue();
            }
        }

        return new ConsensusMap(featureIds, sampleNames, matrix);
    }

    /**
     * Median normalization of intensity matrix.
     */
    public static double[][] medianNormalization(double[][] matrix) {
        int nFeatures = matrix.length;
        int nSamples = matrix[0].length;
        double[][] result = new double[nFeatures][nSamples];

        double[] medians = new double[nSamples];
        for (int j = 0; j < nSamples; j++) {
            List<Double> nonzero = new ArrayList<>();
            for (int i = 0; i < nFeatures; i++) {
                if (matrix[i][j] > 0) nonzero.add(matrix[i][j]);
            }
            if (nonzero.isEmpty()) {
                medians[j] = 1.0;
            } else {
                Collections.sort(nonzero);
                int mid = nonzero.size() / 2;
                medians[j] = nonzero.size() % 2 == 0
                        ? (nonzero.get(mid - 1) + nonzero.get(mid)) / 2.0
                        : nonzero.get(mid);
            }
        }

        // Target median
        List<Double> validMedians = new ArrayList<>();
        for (double m : medians) {
            if (m > 0) validMedians.add(m);
        }
        Collections.sort(validMedians);
        double target = validMedians.isEmpty() ? 1.0 : validMedians.get(validMedians.size() / 2);

        for (int j = 0; j < nSamples; j++) {
            double scale = medians[j] > 0 ? target / medians[j] : 1.0;
            for (int i = 0; i < nFeatures; i++) {
                result[i][j] = matrix[i][j] * scale;
            }
        }

        return result;
    }
}
