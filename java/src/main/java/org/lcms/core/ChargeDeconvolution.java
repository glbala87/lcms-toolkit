package org.lcms.core;

import java.util.*;

/**
 * Charge state deconvolution and isotope pattern detection for mass spectra.
 */
public class ChargeDeconvolution {

    private int minCharge = 1;
    private int maxCharge = 6;
    private double mzTolerance = 0.01;
    private int minPeaks = 2;
    private double minScore = 0.5;
    private double minIntensity = 0.0;
    private double massTolerance = 0.5;

    public ChargeDeconvolution() {}

    public ChargeDeconvolution(int minCharge, int maxCharge, double mzTolerance,
                                int minPeaks, double minScore) {
        this.minCharge = minCharge;
        this.maxCharge = maxCharge;
        this.mzTolerance = mzTolerance;
        this.minPeaks = minPeaks;
        this.minScore = minScore;
    }

    // Setters
    public void setMinCharge(int v) { this.minCharge = v; }
    public void setMaxCharge(int v) { this.maxCharge = v; }
    public void setMzTolerance(double v) { this.mzTolerance = v; }
    public void setMinPeaks(int v) { this.minPeaks = v; }
    public void setMinScore(double v) { this.minScore = v; }
    public void setMinIntensity(double v) { this.minIntensity = v; }
    public void setMassTolerance(double v) { this.massTolerance = v; }

    /**
     * Detect isotope patterns in a spectrum.
     */
    public List<IsotopePattern> detectIsotopePatterns(Spectrum spectrum) {
        double[] mz = spectrum.getMz();
        double[] intensity = spectrum.getIntensity();
        int n = mz.length;

        if (n < minPeaks) return Collections.emptyList();

        // Sort by m/z
        Integer[] sortIdx = new Integer[n];
        for (int i = 0; i < n; i++) sortIdx[i] = i;
        Arrays.sort(sortIdx, (a, b) -> Double.compare(mz[a], mz[b]));

        double[] sortedMz = new double[n];
        double[] sortedInt = new double[n];
        for (int i = 0; i < n; i++) {
            sortedMz[i] = mz[sortIdx[i]];
            sortedInt[i] = intensity[sortIdx[i]];
        }

        // Intensity-sorted indices for seed selection
        Integer[] intOrder = new Integer[n];
        for (int i = 0; i < n; i++) intOrder[i] = i;
        Arrays.sort(intOrder, (a, b) -> Double.compare(sortedInt[b], sortedInt[a]));

        boolean[] used = new boolean[n];
        List<IsotopePattern> patterns = new ArrayList<>();

        for (int seedIdx : intOrder) {
            if (used[seedIdx]) continue;
            if (sortedInt[seedIdx] < minIntensity) continue;

            double seedMz = sortedMz[seedIdx];
            double seedInt = sortedInt[seedIdx];

            IsotopePattern bestPattern = null;
            double bestScore = -1.0;
            List<Integer> bestIndices = null;

            for (int charge = minCharge; charge <= maxCharge; charge++) {
                double expectedSpacing = IsotopePattern.NEUTRON_MASS / charge;

                List<double[]> patternPeaks = new ArrayList<>();
                patternPeaks.add(new double[]{seedMz, seedInt});
                List<Integer> patternIndices = new ArrayList<>();
                patternIndices.add(seedIdx);

                // Forward isotope peaks
                for (int iso = 1; iso <= 10; iso++) {
                    double expectedMz = seedMz + iso * expectedSpacing;
                    int bestCand = -1;
                    double bestDiff = mzTolerance + 1.0;

                    for (int i = 0; i < n; i++) {
                        if (used[i]) continue;
                        double diff = Math.abs(sortedMz[i] - expectedMz);
                        if (diff <= mzTolerance && diff < bestDiff) {
                            bestDiff = diff;
                            bestCand = i;
                        }
                    }

                    if (bestCand < 0) break;
                    patternPeaks.add(new double[]{sortedMz[bestCand], sortedInt[bestCand]});
                    patternIndices.add(bestCand);
                }

                if (patternPeaks.size() < minPeaks) continue;

                // Score against theoretical
                double monoMz = patternPeaks.get(0)[0];
                double neutralMass = (monoMz - IsotopePattern.PROTON_MASS) * charge;
                double[] theoretical = IsotopePattern.averagineDistribution(
                        neutralMass, patternPeaks.size());

                double[] observed = new double[patternPeaks.size()];
                double maxObs = 0;
                for (int i = 0; i < observed.length; i++) {
                    observed[i] = patternPeaks.get(i)[1];
                    if (observed[i] > maxObs) maxObs = observed[i];
                }
                if (maxObs > 0) {
                    for (int i = 0; i < observed.length; i++) observed[i] /= maxObs;
                }

                double score = IsotopePattern.cosineSimilarity(theoretical, observed);

                if (score > bestScore && score >= minScore) {
                    bestScore = score;
                    bestPattern = new IsotopePattern(monoMz, charge, patternPeaks, score);
                    bestPattern.setRetentionTime(spectrum.getRetentionTime());
                    bestIndices = new ArrayList<>(patternIndices);
                }
            }

            if (bestPattern != null) {
                patterns.add(bestPattern);
                for (int idx : bestIndices) used[idx] = true;
            }
        }

        patterns.sort((a, b) -> Double.compare(b.getScore(), a.getScore()));
        return patterns;
    }

    /**
     * Deconvolute spectrum to neutral masses.
     */
    public List<DeconvolutedMass> deconvolute(Spectrum spectrum) {
        List<IsotopePattern> patterns = detectIsotopePatterns(spectrum);
        if (patterns.isEmpty()) return Collections.emptyList();

        List<DeconvolutedMass> masses = new ArrayList<>();

        for (IsotopePattern pattern : patterns) {
            double neutralMass = pattern.getNeutralMass();
            boolean merged = false;

            for (DeconvolutedMass dm : masses) {
                if (Math.abs(dm.getNeutralMass() - neutralMass) <= massTolerance) {
                    dm.merge(pattern);
                    merged = true;
                    break;
                }
            }

            if (!merged) {
                masses.add(new DeconvolutedMass(neutralMass, pattern));
            }
        }

        // Update quality scores
        for (DeconvolutedMass dm : masses) {
            dm.updateQualityScore();
        }

        masses.sort((a, b) -> Double.compare(b.getIntensity(), a.getIntensity()));
        return masses;
    }

    /**
     * Represents a deconvoluted neutral mass.
     */
    public static class DeconvolutedMass {
        private double neutralMass;
        private double intensity;
        private List<Integer> chargeStates;
        private List<IsotopePattern> patterns;
        private double qualityScore;

        public DeconvolutedMass(double neutralMass, IsotopePattern pattern) {
            this.neutralMass = neutralMass;
            this.intensity = pattern.getTotalIntensity();
            this.chargeStates = new ArrayList<>();
            this.chargeStates.add(pattern.getCharge());
            this.patterns = new ArrayList<>();
            this.patterns.add(pattern);
            this.qualityScore = pattern.getScore();
        }

        public void merge(IsotopePattern pattern) {
            intensity += pattern.getTotalIntensity();
            if (!chargeStates.contains(pattern.getCharge())) {
                chargeStates.add(pattern.getCharge());
            }
            patterns.add(pattern);

            // Weighted average mass
            double totalInt = 0, weightedMass = 0;
            for (IsotopePattern p : patterns) {
                double ti = p.getTotalIntensity();
                totalInt += ti;
                weightedMass += p.getNeutralMass() * ti;
            }
            if (totalInt > 0) neutralMass = weightedMass / totalInt;
        }

        public void updateQualityScore() {
            double avgScore = 0;
            for (IsotopePattern p : patterns) avgScore += p.getScore();
            avgScore /= patterns.size();
            double chargeBonus = Math.min(0.2, 0.1 * (chargeStates.size() - 1));
            qualityScore = Math.min(1.0, avgScore + chargeBonus);
        }

        public double getNeutralMass() { return neutralMass; }
        public double getIntensity() { return intensity; }
        public List<Integer> getChargeStates() { return chargeStates; }
        public List<IsotopePattern> getPatterns() { return patterns; }
        public double getQualityScore() { return qualityScore; }

        @Override
        public String toString() {
            return String.format("DeconvolutedMass(mass=%.4f, charges=%s, score=%.3f)",
                    neutralMass, chargeStates, qualityScore);
        }
    }
}
