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
 * Generate extracted ion chromatograms (XICs).
 */
@Command(
    name = "xic",
    description = "Generate extracted ion chromatograms"
)
public class XicCommand implements Callable<Integer> {

    @Parameters(index = "0", description = "Input file (mzML or mzXML)")
    private File inputFile;

    @Parameters(index = "1", description = "Output file (TSV or CSV)")
    private File outputFile;

    @Option(names = {"-m", "--mz"}, required = true, description = "Target m/z value(s)", split = ",")
    private List<Double> targetMzList;

    @Option(names = {"-t", "--tolerance"}, description = "m/z tolerance (default: 0.01)", defaultValue = "0.01")
    private double tolerance;

    @Option(names = {"--ppm"}, description = "Tolerance is in ppm instead of Da")
    private boolean usePpm;

    @Option(names = {"--ms-level"}, description = "MS level to use (default: 1)", defaultValue = "1")
    private int msLevel;

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
        spectra.sort(Comparator.comparingDouble(Spectrum::getRetentionTime));

        System.out.printf("Processing %d MS%d spectra%n", spectra.size(), msLevel);
        System.out.printf("Generating XICs for %d target m/z values%n", targetMzList.size());

        // Generate XICs
        Map<Double, Chromatogram> xics = new LinkedHashMap<>();
        for (Double targetMz : targetMzList) {
            Chromatogram xic = generateXic(spectra, targetMz);
            xics.put(targetMz, xic);
            System.out.printf("  m/z %.4f: apex RT = %.2f s, max intensity = %.2e%n",
                targetMz, xic.getApexRt(), xic.getMaxIntensity());
        }

        // Write output
        String outputName = outputFile.getName().toLowerCase();
        String sep = outputName.endsWith(".csv") ? "," : "\t";

        try (PrintWriter writer = new PrintWriter(new FileWriter(outputFile))) {
            // Header
            List<String> header = new ArrayList<>();
            header.add("rt_seconds");
            for (Double mz : targetMzList) {
                header.add(String.format("xic_%.4f", mz));
            }
            writer.println(String.join(sep, header));

            // Data rows
            // Use first XIC for RT values
            Chromatogram firstXic = xics.values().iterator().next();
            for (int i = 0; i < firstXic.size(); i++) {
                List<String> row = new ArrayList<>();
                row.add(String.format("%.4f", firstXic.getRtAt(i)));

                for (Double mz : targetMzList) {
                    Chromatogram xic = xics.get(mz);
                    row.add(String.format("%.2f", xic.getIntensityAt(i)));
                }
                writer.println(String.join(sep, row));
            }
        }

        System.out.println("Results written to: " + outputFile);
        return 0;
    }

    private Chromatogram generateXic(List<Spectrum> spectra, double targetMz) {
        double[] rt = new double[spectra.size()];
        double[] intensity = new double[spectra.size()];

        for (int i = 0; i < spectra.size(); i++) {
            Spectrum spectrum = spectra.get(i);
            rt[i] = spectrum.getRetentionTime();

            double tol = usePpm ? targetMz * tolerance * 1e-6 : tolerance;
            double mzLow = targetMz - tol;
            double mzHigh = targetMz + tol;

            double sum = 0;
            double[] mzArray = spectrum.getMz();
            double[] intArray = spectrum.getIntensity();

            for (int j = 0; j < mzArray.length; j++) {
                if (mzArray[j] >= mzLow && mzArray[j] <= mzHigh) {
                    sum += intArray[j];
                }
            }
            intensity[i] = sum;
        }

        Chromatogram xic = new Chromatogram(rt, intensity);
        xic.setType(ChromatogramType.XIC);
        xic.setTargetMz(targetMz);
        xic.setMzTolerance(usePpm ? targetMz * tolerance * 1e-6 : tolerance);
        return xic;
    }
}
