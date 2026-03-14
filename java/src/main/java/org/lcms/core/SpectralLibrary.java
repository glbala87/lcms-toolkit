package org.lcms.core;

import java.io.*;
import java.util.*;

/**
 * Spectral library for MS/MS matching and searching.
 */
public class SpectralLibrary {

    private final List<LibrarySpectrum> entries = new ArrayList<>();

    public int size() { return entries.size(); }

    public LibrarySpectrum get(int index) { return entries.get(index); }

    public void addSpectrum(String name, double precursorMz,
                            double[] mz, double[] intensity,
                            Map<String, String> metadata) {
        entries.add(new LibrarySpectrum(name, precursorMz, mz, intensity,
                metadata != null ? metadata : new HashMap<>()));
    }

    /**
     * Load spectra from an MGF file.
     */
    public int loadMGF(String filepath) throws IOException {
        int count = 0;
        try (BufferedReader reader = new BufferedReader(new FileReader(filepath))) {
            boolean inIons = false;
            String name = "";
            double precursorMz = 0;
            Map<String, String> metadata = new HashMap<>();
            List<Double> mzList = new ArrayList<>();
            List<Double> intList = new ArrayList<>();

            String line;
            while ((line = reader.readLine()) != null) {
                line = line.trim();

                if ("BEGIN IONS".equals(line)) {
                    inIons = true;
                    name = "";
                    precursorMz = 0;
                    metadata = new HashMap<>();
                    mzList = new ArrayList<>();
                    intList = new ArrayList<>();
                } else if ("END IONS".equals(line)) {
                    if (!mzList.isEmpty()) {
                        addSpectrum(name, precursorMz,
                                toArray(mzList), toArray(intList), metadata);
                        count++;
                    }
                    inIons = false;
                } else if (inIons) {
                    int eqIdx = line.indexOf('=');
                    if (eqIdx >= 0) {
                        String key = line.substring(0, eqIdx).trim().toUpperCase();
                        String value = line.substring(eqIdx + 1).trim();
                        switch (key) {
                            case "TITLE": name = value; break;
                            case "PEPMASS":
                                precursorMz = Double.parseDouble(value.split("\\s+")[0]);
                                break;
                            default: metadata.put(key.toLowerCase(), value);
                        }
                    } else {
                        String[] parts = line.split("\\s+");
                        if (parts.length >= 2) {
                            try {
                                mzList.add(Double.parseDouble(parts[0]));
                                intList.add(Double.parseDouble(parts[1]));
                            } catch (NumberFormatException ignored) {}
                        }
                    }
                }
            }
        }
        return count;
    }

    /**
     * Load spectra from an MSP file.
     */
    public int loadMSP(String filepath) throws IOException {
        int count = 0;
        try (BufferedReader reader = new BufferedReader(new FileReader(filepath))) {
            String name = "";
            double precursorMz = 0;
            Map<String, String> metadata = new HashMap<>();
            List<Double> mzList = new ArrayList<>();
            List<Double> intList = new ArrayList<>();
            boolean readingPeaks = false;

            String line;
            while ((line = reader.readLine()) != null) {
                line = line.trim();

                if (line.isEmpty()) {
                    if (!mzList.isEmpty()) {
                        addSpectrum(name, precursorMz,
                                toArray(mzList), toArray(intList), metadata);
                        count++;
                    }
                    name = "";
                    precursorMz = 0;
                    metadata = new HashMap<>();
                    mzList = new ArrayList<>();
                    intList = new ArrayList<>();
                    readingPeaks = false;
                    continue;
                }

                if (!readingPeaks) {
                    int colonIdx = line.indexOf(':');
                    if (colonIdx >= 0) {
                        String key = line.substring(0, colonIdx).trim();
                        String value = line.substring(colonIdx + 1).trim();
                        String keyUpper = key.toUpperCase();

                        if ("NAME".equals(keyUpper)) {
                            name = value;
                        } else if (keyUpper.contains("PRECURSORMZ")) {
                            precursorMz = Double.parseDouble(value);
                        } else if (keyUpper.contains("NUM") && keyUpper.contains("PEAK")) {
                            readingPeaks = true;
                        } else {
                            metadata.put(key.toLowerCase(), value);
                        }
                    }
                } else {
                    String[] parts = line.split("[\\t ]+");
                    if (parts.length >= 2) {
                        try {
                            mzList.add(Double.parseDouble(parts[0]));
                            intList.add(Double.parseDouble(parts[1]));
                        } catch (NumberFormatException ignored) {}
                    }
                }
            }

            if (!mzList.isEmpty()) {
                addSpectrum(name, precursorMz,
                        toArray(mzList), toArray(intList), metadata);
                count++;
            }
        }
        return count;
    }

    /**
     * Save library to MGF format.
     */
    public void saveMGF(String filepath) throws IOException {
        try (PrintWriter writer = new PrintWriter(new FileWriter(filepath))) {
            for (LibrarySpectrum entry : entries) {
                writer.println("BEGIN IONS");
                writer.printf("TITLE=%s%n", entry.name);
                writer.printf("PEPMASS=%.6f%n", entry.precursorMz);
                for (Map.Entry<String, String> meta : entry.metadata.entrySet()) {
                    writer.printf("%s=%s%n", meta.getKey().toUpperCase(), meta.getValue());
                }
                for (int i = 0; i < entry.mz.length; i++) {
                    writer.printf("%.6f %.4f%n", entry.mz[i], entry.intensity[i]);
                }
                writer.println("END IONS");
                writer.println();
            }
        }
    }

    /**
     * Search library for matching spectra.
     */
    public List<SpectralMatch> search(double[] queryMz, double[] queryIntensity,
                                       double queryPrecursorMz, String method,
                                       double tolerance, double minScore, int topN,
                                       double precursorTolerance) {
        List<SpectralMatch> results = new ArrayList<>();

        for (LibrarySpectrum entry : entries) {
            if (precursorTolerance > 0 && queryPrecursorMz > 0) {
                if (Math.abs(entry.precursorMz - queryPrecursorMz) > precursorTolerance) {
                    continue;
                }
            }

            double[] scoreAndMatched = computeSimilarity(
                    queryMz, queryIntensity, queryPrecursorMz,
                    entry.mz, entry.intensity, entry.precursorMz,
                    method, tolerance);

            if (scoreAndMatched[0] >= minScore) {
                results.add(new SpectralMatch(
                        entry.name, scoreAndMatched[0],
                        (int) scoreAndMatched[1], entry.precursorMz, entry.metadata));
            }
        }

        results.sort((a, b) -> Double.compare(b.getScore(), a.getScore()));
        if (results.size() > topN) {
            results = results.subList(0, topN);
        }

        return results;
    }

    private double[] computeSimilarity(double[] mz1, double[] int1, double pre1,
                                        double[] mz2, double[] int2, double pre2,
                                        String method, double tolerance) {
        switch (method) {
            case "modified_cosine":
                return modifiedCosineSimilarity(mz1, int1, pre1, mz2, int2, pre2, tolerance);
            case "cosine":
            default:
                return cosineSimilarity(mz1, int1, mz2, int2, tolerance);
        }
    }

    public static double[] cosineSimilarity(double[] mz1, double[] int1,
                                             double[] mz2, double[] int2,
                                             double tolerance) {
        List<double[]> matched = matchPeaks(mz1, int1, mz2, int2, tolerance);
        if (matched.isEmpty()) return new double[]{0, 0};

        double dot = 0;
        for (double[] pair : matched) dot += pair[0] * pair[1];

        double norm1 = vectorNorm(int1);
        double norm2 = vectorNorm(int2);

        if (norm1 == 0 || norm2 == 0) return new double[]{0, 0};

        double score = Math.min(dot / (norm1 * norm2), 1.0);
        return new double[]{score, matched.size()};
    }

    public static double[] modifiedCosineSimilarity(
            double[] mz1, double[] int1, double pre1,
            double[] mz2, double[] int2, double pre2, double tolerance) {

        double massDiff = pre1 - pre2;
        List<double[]> direct = matchPeaks(mz1, int1, mz2, int2, tolerance);

        double[] shiftedMz2 = new double[mz2.length];
        for (int i = 0; i < mz2.length; i++) shiftedMz2[i] = mz2[i] + massDiff;
        List<double[]> shifted = matchPeaks(mz1, int1, shiftedMz2, int2, tolerance);

        List<double[]> all = new ArrayList<>(direct);
        all.addAll(shifted);
        if (all.isEmpty()) return new double[]{0, 0};

        double dot = 0;
        for (double[] pair : all) dot += pair[0] * pair[1];

        double norm1 = vectorNorm(int1);
        double norm2 = vectorNorm(int2);

        if (norm1 == 0 || norm2 == 0) return new double[]{0, 0};

        double score = Math.min(dot / (norm1 * norm2), 1.0);
        return new double[]{score, all.size()};
    }

    private static List<double[]> matchPeaks(double[] mz1, double[] int1,
                                              double[] mz2, double[] int2,
                                              double tolerance) {
        List<double[]> matched = new ArrayList<>();
        boolean[] used2 = new boolean[mz2.length];

        for (int i = 0; i < mz1.length; i++) {
            int bestJ = -1;
            double bestDiff = tolerance + 1;

            for (int j = 0; j < mz2.length; j++) {
                if (used2[j]) continue;
                double diff = Math.abs(mz1[i] - mz2[j]);
                if (diff <= tolerance && diff < bestDiff) {
                    bestDiff = diff;
                    bestJ = j;
                }
            }

            if (bestJ >= 0) {
                matched.add(new double[]{int1[i], int2[bestJ]});
                used2[bestJ] = true;
            }
        }

        return matched;
    }

    private static double vectorNorm(double[] v) {
        double sum = 0;
        for (double x : v) sum += x * x;
        return Math.sqrt(sum);
    }

    private static double[] toArray(List<Double> list) {
        double[] arr = new double[list.size()];
        for (int i = 0; i < list.size(); i++) arr[i] = list.get(i);
        return arr;
    }

    /**
     * Library spectrum entry.
     */
    public static class LibrarySpectrum {
        public final String name;
        public final double precursorMz;
        public final double[] mz;
        public final double[] intensity;
        public final Map<String, String> metadata;

        public LibrarySpectrum(String name, double precursorMz,
                                double[] mz, double[] intensity,
                                Map<String, String> metadata) {
            this.name = name;
            this.precursorMz = precursorMz;
            this.mz = mz;
            this.intensity = intensity;
            this.metadata = metadata;
        }

        public int getNumPeaks() { return mz.length; }
    }
}
