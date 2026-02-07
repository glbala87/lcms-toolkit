package org.lcms.tools;

import org.lcms.core.*;
import org.lcms.io.MzMLReader;
import org.lcms.io.MzXMLReader;
import picocli.CommandLine.Command;
import picocli.CommandLine.Parameters;
import picocli.CommandLine.Option;

import java.io.File;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.Callable;

/**
 * Display information about an LC-MS file.
 */
@Command(
    name = "info",
    description = "Display information about an LC-MS data file"
)
public class InfoCommand implements Callable<Integer> {

    @Parameters(index = "0", description = "Input file (mzML or mzXML)")
    private File inputFile;

    @Option(names = {"-d", "--detailed"}, description = "Show detailed statistics")
    private boolean detailed;

    @Override
    public Integer call() throws Exception {
        if (!inputFile.exists()) {
            System.err.println("Error: File not found: " + inputFile);
            return 1;
        }

        String filename = inputFile.getName().toLowerCase();
        MSExperiment experiment;

        try {
            if (filename.endsWith(".mzml")) {
                MzMLReader reader = new MzMLReader();
                experiment = reader.read(inputFile.getAbsolutePath());
            } else if (filename.endsWith(".mzxml")) {
                MzXMLReader reader = new MzXMLReader();
                experiment = reader.read(inputFile.getAbsolutePath());
            } else {
                System.err.println("Error: Unsupported file format. Use .mzML or .mzXML");
                return 1;
            }
        } catch (Exception e) {
            System.err.println("Error reading file: " + e.getMessage());
            return 1;
        }

        // Print basic info
        System.out.println("=== File Information ===");
        System.out.println("File: " + inputFile.getAbsolutePath());
        System.out.println("Size: " + formatBytes(inputFile.length()));

        if (experiment.getInstrumentModel() != null) {
            System.out.println("Instrument: " + experiment.getInstrumentModel());
        }
        if (experiment.getDateTime() != null) {
            System.out.println("Date: " + experiment.getDateTime());
        }

        System.out.println();
        System.out.println("=== Data Summary ===");
        System.out.println("Spectra: " + experiment.getSpectrumCount());
        System.out.println("Chromatograms: " + experiment.getChromatogramCount());
        System.out.println("Total data points: " + experiment.getTotalDataPoints());
        System.out.printf("Average spectrum size: %.1f points%n", experiment.getAverageSpectrumSize());

        // MS levels breakdown
        Map<Integer, Integer> msLevels = new HashMap<>();
        for (Spectrum s : experiment.getSpectra()) {
            msLevels.merge(s.getMsLevel(), 1, Integer::sum);
        }

        System.out.println();
        System.out.println("=== MS Levels ===");
        for (Map.Entry<Integer, Integer> entry : msLevels.entrySet()) {
            System.out.printf("MS%d: %d spectra%n", entry.getKey(), entry.getValue());
        }

        // Ranges
        double[] mzRange = experiment.getMzRange();
        double[] rtRange = experiment.getRtRange();

        System.out.println();
        System.out.println("=== Ranges ===");
        System.out.printf("m/z: %.4f - %.4f%n", mzRange[0], mzRange[1]);
        System.out.printf("RT: %.2f - %.2f seconds%n", rtRange[0], rtRange[1]);
        System.out.printf("RT: %.2f - %.2f minutes%n", rtRange[0] / 60, rtRange[1] / 60);

        // Detailed statistics
        if (detailed && experiment.hasSpectra()) {
            System.out.println();
            System.out.println("=== Detailed Statistics ===");

            // TIC statistics
            double maxTic = 0;
            double minTic = Double.POSITIVE_INFINITY;
            double sumTic = 0;

            for (Spectrum s : experiment.getSpectra()) {
                if (s.getMsLevel() == 1) {
                    if (s.getTic() > maxTic) maxTic = s.getTic();
                    if (s.getTic() < minTic) minTic = s.getTic();
                    sumTic += s.getTic();
                }
            }

            int ms1Count = experiment.countSpectraByLevel(1);
            if (ms1Count > 0) {
                System.out.printf("MS1 TIC range: %.2e - %.2e%n", minTic, maxTic);
                System.out.printf("MS1 TIC mean: %.2e%n", sumTic / ms1Count);
            }

            // Spectrum size statistics
            int minSize = Integer.MAX_VALUE;
            int maxSize = 0;

            for (Spectrum s : experiment.getSpectra()) {
                if (s.size() < minSize) minSize = s.size();
                if (s.size() > maxSize) maxSize = s.size();
            }

            System.out.printf("Spectrum size range: %d - %d points%n", minSize, maxSize);

            // Chromatogram info
            if (experiment.hasChromatograms()) {
                System.out.println();
                System.out.println("=== Chromatograms ===");
                for (Chromatogram c : experiment.getChromatograms()) {
                    System.out.printf("  %s: %s (%d points)%n",
                        c.getType(), c.getNativeId(), c.size());
                }
            }
        }

        return 0;
    }

    private String formatBytes(long bytes) {
        if (bytes < 1024) return bytes + " B";
        if (bytes < 1024 * 1024) return String.format("%.1f KB", bytes / 1024.0);
        if (bytes < 1024 * 1024 * 1024) return String.format("%.1f MB", bytes / (1024.0 * 1024));
        return String.format("%.1f GB", bytes / (1024.0 * 1024 * 1024));
    }
}
