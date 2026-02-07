package org.lcms.tools;

import org.lcms.core.*;
import org.lcms.io.MzMLReader;
import org.lcms.io.MzXMLReader;
import picocli.CommandLine.Command;
import picocli.CommandLine.Parameters;
import picocli.CommandLine.Option;

import java.io.*;
import java.util.concurrent.Callable;

/**
 * Convert LC-MS data files between formats.
 */
@Command(
    name = "convert",
    description = "Convert LC-MS data files between formats"
)
public class ConvertCommand implements Callable<Integer> {

    @Parameters(index = "0", description = "Input file")
    private File inputFile;

    @Parameters(index = "1", description = "Output file")
    private File outputFile;

    @Option(names = {"--ms-level"}, description = "Only convert spectra at this MS level")
    private Integer msLevel;

    @Option(names = {"--rt-min"}, description = "Minimum retention time (seconds)")
    private Double rtMin;

    @Option(names = {"--rt-max"}, description = "Maximum retention time (seconds)")
    private Double rtMax;

    @Option(names = {"--max-spectra"}, description = "Maximum number of spectra to convert")
    private Integer maxSpectra;

    @Override
    public Integer call() throws Exception {
        if (!inputFile.exists()) {
            System.err.println("Error: Input file not found: " + inputFile);
            return 1;
        }

        String inputName = inputFile.getName().toLowerCase();
        String outputName = outputFile.getName().toLowerCase();

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

        System.out.printf("Loaded %d spectra, %d chromatograms%n",
            experiment.getSpectrumCount(), experiment.getChromatogramCount());

        // Apply filters
        if (msLevel != null || rtMin != null || rtMax != null) {
            MSExperiment filtered = new MSExperiment();
            filtered.setSourceFile(experiment.getSourceFile());
            filtered.setInstrumentModel(experiment.getInstrumentModel());

            int count = 0;
            for (Spectrum s : experiment.getSpectra()) {
                if (maxSpectra != null && count >= maxSpectra) break;

                if (msLevel != null && s.getMsLevel() != msLevel) continue;
                if (rtMin != null && s.getRetentionTime() < rtMin) continue;
                if (rtMax != null && s.getRetentionTime() > rtMax) continue;

                filtered.addSpectrum(s);
                count++;
            }

            for (Chromatogram c : experiment.getChromatograms()) {
                filtered.addChromatogram(c);
            }

            experiment = filtered;
            System.out.printf("After filtering: %d spectra%n", experiment.getSpectrumCount());
        }

        // Write output
        System.out.println("Writing: " + outputFile);
        try {
            if (outputName.endsWith(".tsv") || outputName.endsWith(".txt")) {
                writeTsv(experiment, outputFile);
            } else if (outputName.endsWith(".csv")) {
                writeCsv(experiment, outputFile);
            } else {
                System.err.println("Error: Unsupported output format. Use .tsv, .txt, or .csv");
                return 1;
            }
        } catch (Exception e) {
            System.err.println("Error writing file: " + e.getMessage());
            return 1;
        }

        System.out.println("Done!");
        return 0;
    }

    private void writeTsv(MSExperiment exp, File output) throws IOException {
        try (PrintWriter writer = new PrintWriter(new FileWriter(output))) {
            // Write summary header
            writer.println("# LCMS Export");
            writer.println("# Source: " + exp.getSourceFile());
            writer.println("# Spectra: " + exp.getSpectrumCount());
            writer.println();

            // Write spectrum summary
            writer.println("spectrum_index\tms_level\trt_seconds\ttic\tbase_peak_mz\tbase_peak_intensity\tnum_peaks");
            for (Spectrum s : exp.getSpectra()) {
                writer.printf("%d\t%d\t%.4f\t%.2f\t%.6f\t%.2f\t%d%n",
                    s.getIndex(),
                    s.getMsLevel(),
                    s.getRetentionTime(),
                    s.getTic(),
                    s.getBasePeakMz(),
                    s.getBasePeakIntensity(),
                    s.size()
                );
            }
        }
    }

    private void writeCsv(MSExperiment exp, File output) throws IOException {
        try (PrintWriter writer = new PrintWriter(new FileWriter(output))) {
            writer.println("spectrum_index,ms_level,rt_seconds,tic,base_peak_mz,base_peak_intensity,num_peaks");
            for (Spectrum s : exp.getSpectra()) {
                writer.printf("%d,%d,%.4f,%.2f,%.6f,%.2f,%d%n",
                    s.getIndex(),
                    s.getMsLevel(),
                    s.getRetentionTime(),
                    s.getTic(),
                    s.getBasePeakMz(),
                    s.getBasePeakIntensity(),
                    s.size()
                );
            }
        }
    }
}
