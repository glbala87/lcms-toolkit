package org.lcms.core;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Container for a complete LC-MS experiment.
 */
public class MSExperiment {
    private List<Spectrum> spectra;
    private List<Chromatogram> chromatograms;
    private String sourceFile;
    private String dateTime;
    private String instrumentModel;
    private String instrumentSerial;
    private String software;
    private Map<String, String> metadata;

    public MSExperiment() {
        this.spectra = new ArrayList<>();
        this.chromatograms = new ArrayList<>();
        this.metadata = new HashMap<>();
    }

    // Spectra access
    public int getSpectrumCount() { return spectra.size(); }
    public boolean hasSpectra() { return !spectra.isEmpty(); }

    public Spectrum getSpectrum(int index) { return spectra.get(index); }
    public List<Spectrum> getSpectra() { return spectra; }

    public void addSpectrum(Spectrum spectrum) {
        spectrum.setIndex(spectra.size());
        spectra.add(spectrum);
    }

    public void reserveSpectra(int n) {
        ((ArrayList<Spectrum>) spectra).ensureCapacity(n);
    }

    // Chromatograms access
    public int getChromatogramCount() { return chromatograms.size(); }
    public boolean hasChromatograms() { return !chromatograms.isEmpty(); }

    public Chromatogram getChromatogram(int index) { return chromatograms.get(index); }
    public List<Chromatogram> getChromatograms() { return chromatograms; }

    public void addChromatogram(Chromatogram chrom) {
        chrom.setIndex(chromatograms.size());
        chromatograms.add(chrom);
    }

    // Filtering
    public List<Spectrum> getSpectraByLevel(int level) {
        return spectra.stream()
            .filter(s -> s.getMsLevel() == level)
            .collect(Collectors.toList());
    }

    public int countSpectraByLevel(int level) {
        return (int) spectra.stream()
            .filter(s -> s.getMsLevel() == level)
            .count();
    }

    public List<Spectrum> getSpectraInRtRange(double rtMin, double rtMax) {
        return spectra.stream()
            .filter(s -> s.getRetentionTime() >= rtMin && s.getRetentionTime() <= rtMax)
            .collect(Collectors.toList());
    }

    public Optional<Spectrum> findSpectrumByNativeId(String nativeId) {
        return spectra.stream()
            .filter(s -> nativeId.equals(s.getNativeId()))
            .findFirst();
    }

    public Optional<Spectrum> findSpectrumByRt(double rt, int level) {
        return spectra.stream()
            .filter(s -> level == 0 || s.getMsLevel() == level)
            .min(Comparator.comparingDouble(s -> Math.abs(s.getRetentionTime() - rt)));
    }

    // Chromatogram generation
    public Chromatogram generateTIC(int level) {
        List<Spectrum> filtered = level == 0 ? spectra : getSpectraByLevel(level);

        double[] rt = new double[filtered.size()];
        double[] intensity = new double[filtered.size()];

        for (int i = 0; i < filtered.size(); i++) {
            rt[i] = filtered.get(i).getRetentionTime();
            intensity[i] = filtered.get(i).getTic();
        }

        Chromatogram tic = new Chromatogram(rt, intensity);
        tic.setType(ChromatogramType.TIC);
        tic.setNativeId("TIC");
        tic.sortByRt();
        return tic;
    }

    public Chromatogram generateBPC(int level) {
        List<Spectrum> filtered = level == 0 ? spectra : getSpectraByLevel(level);

        double[] rt = new double[filtered.size()];
        double[] intensity = new double[filtered.size()];

        for (int i = 0; i < filtered.size(); i++) {
            rt[i] = filtered.get(i).getRetentionTime();
            intensity[i] = filtered.get(i).getBasePeakIntensity();
        }

        Chromatogram bpc = new Chromatogram(rt, intensity);
        bpc.setType(ChromatogramType.BPC);
        bpc.setNativeId("BPC");
        bpc.sortByRt();
        return bpc;
    }

    // Ranges
    public double[] getMzRange() {
        if (spectra.isEmpty()) return new double[]{0, 0};

        double min = Double.POSITIVE_INFINITY;
        double max = Double.NEGATIVE_INFINITY;

        for (Spectrum s : spectra) {
            if (s.size() > 0) {
                if (s.getMzMin() < min) min = s.getMzMin();
                if (s.getMzMax() > max) max = s.getMzMax();
            }
        }
        return new double[]{min, max};
    }

    public double[] getRtRange() {
        if (spectra.isEmpty()) return new double[]{0, 0};

        double min = Double.POSITIVE_INFINITY;
        double max = Double.NEGATIVE_INFINITY;

        for (Spectrum s : spectra) {
            if (s.getRetentionTime() < min) min = s.getRetentionTime();
            if (s.getRetentionTime() > max) max = s.getRetentionTime();
        }
        return new double[]{min, max};
    }

    // Sorting
    public void sortSpectraByRt() {
        spectra.sort(Comparator.comparingDouble(Spectrum::getRetentionTime));
        for (int i = 0; i < spectra.size(); i++) {
            spectra.get(i).setIndex(i);
        }
    }

    // Clear
    public void clear() {
        spectra.clear();
        chromatograms.clear();
        metadata.clear();
    }

    // Metadata
    public String getSourceFile() { return sourceFile; }
    public void setSourceFile(String sourceFile) { this.sourceFile = sourceFile; }

    public String getDateTime() { return dateTime; }
    public void setDateTime(String dateTime) { this.dateTime = dateTime; }

    public String getInstrumentModel() { return instrumentModel; }
    public void setInstrumentModel(String model) { this.instrumentModel = model; }

    public String getInstrumentSerial() { return instrumentSerial; }
    public void setInstrumentSerial(String serial) { this.instrumentSerial = serial; }

    public String getSoftware() { return software; }
    public void setSoftware(String software) { this.software = software; }

    public Map<String, String> getMetadata() { return metadata; }

    // Statistics
    public long getTotalDataPoints() {
        return spectra.stream().mapToLong(Spectrum::size).sum();
    }

    public double getAverageSpectrumSize() {
        if (spectra.isEmpty()) return 0;
        return (double) getTotalDataPoints() / spectra.size();
    }

    @Override
    public String toString() {
        return String.format("MSExperiment(spectra=%d, chromatograms=%d, source=%s)",
            spectra.size(), chromatograms.size(), sourceFile);
    }
}
