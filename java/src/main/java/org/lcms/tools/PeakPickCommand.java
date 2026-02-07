package org.lcms.tools;

import org.lcms.core.*;
import org.lcms.io.MzMLReader;
import org.lcms.io.MzXMLReader;
import picocli.CommandLine.Command;
import picocli.CommandLine.Parameters;
import picocli.CommandLine.Option;

import java.io.*;
import java.util.*;
import java.util.concurrent.Callable;

/**
 * Pick peaks from LC-MS spectra.
 */
@Command(
    name = "peaks",
    description = "Pick peaks from spectra and export results"
)
public class PeakPickCommand implements Callable<Integer> {

    @Parameters(index = "0", description = "Input file (mzML or mzXML)")
    private File inputFile;

    @Parameters(index = "1", description = "Output file (TSV or CSV)")
    private File outputFile;

    @Option(names = {"--ms-level"}, description = "MS level to process (default: 1)", defaultValue = "1")
    private int msLevel;

    @Option(names = {"--min-snr"}, description = "Minimum signal-to-noise ratio (default: 3.0)", defaultValue = "3.0")
    private double minSnr;

    @Option(names = {"--min-intensity"}, description = "Minimum absolute intensity (default: 0)")
    private double minIntensity;

    @Option(names = {"--window-size"}, description = "Peak detection window size (default: 3)", defaultValue = "3")
    private int windowSize;

    @Option(names = {"--top-n"}, description = "Only output top N peaks per spectrum (0 = all)")
    private int topN;

    @Override
    public Integer call() throws Exception {
        if (!inputFile.exists()) {
            System.err.println("Error: Input file not found: " + inputFile);
            return 1;
        }

        String inputName = inputFile.getName().toLowerCase();
        MSExperiment experiment;

        // Read input
        System.out.println("Reading: " + inputFile);
        try {
            if (inputName.endsWith(".mzml")) {
                MzMLReader reader = new MzMLReader();
                experiment = reader.read(inputFile.getAbsolutePath());
            } else if (inputName.endsWith(".mzxml")) {
                MzXMLReader reader = new MzXMLReader();
                experiment = reader.read(inputFile.getAbsolutePath());
            } else {
                System.err.println("Error: Unsupported input format");
                return 1;
            }
        } catch (Exception e) {
            System.err.println("Error reading file: " + e.getMessage());
            return 1;
        }

        List<Spectrum> spectra = experiment.getSpectraByLevel(msLevel);
        System.out.printf("Processing %d MS%d spectra%n", spectra.size(), msLevel);

        // Peak picking
        List<Peak> allPeaks = new ArrayList<>();
        int totalPeaks = 0;

        for (Spectrum spectrum : spectra) {
            List<Peak> peaks = pickPeaks(spectrum);
            totalPeaks += peaks.size();

            // Apply top-N filter
            if (topN > 0 && peaks.size() > topN) {
                peaks.sort(Comparator.comparingDouble(Peak::getIntensity).reversed());
                peaks = peaks.subList(0, topN);
            }

            allPeaks.addAll(peaks);
        }

        System.out.printf("Found %d peaks total%n", totalPeaks);
        System.out.printf("Output %d peaks%n", allPeaks.size());

        // Write output
        String outputName = outputFile.getName().toLowerCase();
        try (PrintWriter writer = new PrintWriter(new FileWriter(outputFile))) {
            String sep = outputName.endsWith(".csv") ? "," : "\t";

            writer.println(String.join(sep, Arrays.asList(
                "spectrum_index", "rt_seconds", "mz", "intensity", "snr", "area", "fwhm"
            )));

            for (Peak peak : allPeaks) {
                writer.printf("%d%s%.4f%s%.6f%s%.2f%s%.2f%s%.2f%s%.6f%n",
                    peak.getSpectrumIndex(), sep,
                    peak.getRt(), sep,
                    peak.getMz(), sep,
                    peak.getIntensity(), sep,
                    peak.getSnr(), sep,
                    peak.getArea(), sep,
                    peak.getFwhmMz()
                );
            }
        }

        System.out.println("Results written to: " + outputFile);
        return 0;
    }

    private List<Peak> pickPeaks(Spectrum spectrum) {
        List<Peak> peaks = new ArrayList<>();
        double[] mz = spectrum.getMz();
        double[] intensity = spectrum.getIntensity();

        if (mz.length < 3) return peaks;

        // Estimate noise using MAD
        double[] sorted = intensity.clone();
        Arrays.sort(sorted);
        double median = sorted[sorted.length / 2];
        double[] deviations = new double[sorted.length];
        for (int i = 0; i < sorted.length; i++) {
            deviations[i] = Math.abs(sorted[i] - median);
        }
        Arrays.sort(deviations);
        double mad = deviations[deviations.length / 2];
        double noise = mad * 1.4826;
        if (noise <= 0) noise = 1.0;

        int halfWindow = windowSize / 2;

        // Find local maxima
        for (int i = halfWindow; i < mz.length - halfWindow; i++) {
            // Check if local maximum
            boolean isMax = true;
            for (int j = i - halfWindow; j <= i + halfWindow; j++) {
                if (j != i && intensity[j] >= intensity[i]) {
                    isMax = false;
                    break;
                }
            }

            if (!isMax) continue;
            if (intensity[i] < minIntensity) continue;

            double snr = intensity[i] / noise;
            if (snr < minSnr) continue;

            // Find boundaries
            int left = i;
            while (left > 0 && intensity[left - 1] < intensity[left]) {
                left--;
            }

            int right = i;
            while (right < mz.length - 1 && intensity[right + 1] < intensity[right]) {
                right++;
            }

            // Create peak
            Peak peak = new Peak();
            peak.setMz(mz[i]);
            peak.setRt(spectrum.getRetentionTime());
            peak.setIntensity(intensity[i]);
            peak.setSnr(snr);
            peak.setSpectrumIndex(spectrum.getIndex());
            peak.setMzLeft(mz[left]);
            peak.setMzRight(mz[right]);

            // Calculate area (trapezoidal)
            double area = 0;
            for (int j = left; j < right; j++) {
                area += 0.5 * (intensity[j] + intensity[j + 1]) * (mz[j + 1] - mz[j]);
            }
            peak.setArea(area);

            // Calculate FWHM
            double halfMax = intensity[i] / 2;
            double fwhmLeft = mz[i];
            double fwhmRight = mz[i];

            for (int j = i; j >= left; j--) {
                if (intensity[j] <= halfMax) {
                    if (j < i && intensity[j + 1] > halfMax) {
                        double t = (halfMax - intensity[j]) / (intensity[j + 1] - intensity[j]);
                        fwhmLeft = mz[j] + t * (mz[j + 1] - mz[j]);
                    } else {
                        fwhmLeft = mz[j];
                    }
                    break;
                }
            }

            for (int j = i; j <= right; j++) {
                if (intensity[j] <= halfMax) {
                    if (j > i && intensity[j - 1] > halfMax) {
                        double t = (halfMax - intensity[j - 1]) / (intensity[j] - intensity[j - 1]);
                        fwhmRight = mz[j - 1] + t * (mz[j] - mz[j - 1]);
                    } else {
                        fwhmRight = mz[j];
                    }
                    break;
                }
            }

            peak.setFwhmMz(fwhmRight - fwhmLeft);

            peaks.add(peak);
        }

        return peaks;
    }

    private static class Peak {
        private double mz;
        private double rt;
        private double intensity;
        private double snr;
        private double area;
        private double fwhmMz;
        private double mzLeft;
        private double mzRight;
        private int spectrumIndex;

        public double getMz() { return mz; }
        public void setMz(double mz) { this.mz = mz; }

        public double getRt() { return rt; }
        public void setRt(double rt) { this.rt = rt; }

        public double getIntensity() { return intensity; }
        public void setIntensity(double intensity) { this.intensity = intensity; }

        public double getSnr() { return snr; }
        public void setSnr(double snr) { this.snr = snr; }

        public double getArea() { return area; }
        public void setArea(double area) { this.area = area; }

        public double getFwhmMz() { return fwhmMz; }
        public void setFwhmMz(double fwhmMz) { this.fwhmMz = fwhmMz; }

        public double getMzLeft() { return mzLeft; }
        public void setMzLeft(double mzLeft) { this.mzLeft = mzLeft; }

        public double getMzRight() { return mzRight; }
        public void setMzRight(double mzRight) { this.mzRight = mzRight; }

        public int getSpectrumIndex() { return spectrumIndex; }
        public void setSpectrumIndex(int spectrumIndex) { this.spectrumIndex = spectrumIndex; }
    }
}
