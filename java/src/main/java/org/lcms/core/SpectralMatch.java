package org.lcms.core;

import java.util.Map;

/**
 * Result of a spectral library search.
 */
public class SpectralMatch {
    private final String name;
    private final double score;
    private final int matchedPeaks;
    private final double precursorMz;
    private final Map<String, String> metadata;

    public SpectralMatch(String name, double score, int matchedPeaks,
                          double precursorMz, Map<String, String> metadata) {
        this.name = name;
        this.score = score;
        this.matchedPeaks = matchedPeaks;
        this.precursorMz = precursorMz;
        this.metadata = metadata;
    }

    public String getName() { return name; }
    public double getScore() { return score; }
    public int getMatchedPeaks() { return matchedPeaks; }
    public double getPrecursorMz() { return precursorMz; }
    public Map<String, String> getMetadata() { return metadata; }

    @Override
    public String toString() {
        return String.format("SpectralMatch(name='%s', score=%.4f, matched=%d)",
                name, score, matchedPeaks);
    }
}
