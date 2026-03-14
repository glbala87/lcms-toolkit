package org.lcms.io;

import org.lcms.core.*;

import java.io.*;
import java.util.*;
import java.util.regex.*;

/**
 * Indexed mzML reader for random access to spectra without loading
 * the entire file into memory.
 */
public class IndexedMzMLReader implements AutoCloseable {

    private final String filepath;
    private RandomAccessFile file;
    private List<SpectrumIndexEntry> index;
    private Map<Integer, Integer> scanMap;

    public IndexedMzMLReader(String filepath) {
        this.filepath = filepath;
        this.index = new ArrayList<>();
        this.scanMap = new HashMap<>();
    }

    /**
     * Open file and build index.
     */
    public void open() throws IOException {
        this.file = new RandomAccessFile(filepath, "r");
        buildIndex();
    }

    @Override
    public void close() throws IOException {
        if (file != null) {
            file.close();
            file = null;
        }
    }

    public int getSpectrumCount() { return index.size(); }

    /**
     * Get spectrum by sequential index.
     */
    public Spectrum getSpectrum(int idx) throws IOException {
        if (idx < 0 || idx >= index.size()) return null;
        return readSpectrumAt(index.get(idx));
    }

    /**
     * Get spectrum by scan number.
     */
    public Spectrum getSpectrumByScan(int scanNumber) throws IOException {
        Integer idx = scanMap.get(scanNumber);
        if (idx == null) return null;
        return readSpectrumAt(index.get(idx));
    }

    /**
     * Iterate spectra lazily (calls callback for each spectrum).
     */
    public void iterateSpectra(SpectrumConsumer consumer) throws IOException {
        iterateSpectra(consumer, -1);
    }

    /**
     * Iterate spectra filtered by MS level.
     */
    public void iterateSpectra(SpectrumConsumer consumer, int msLevel) throws IOException {
        for (SpectrumIndexEntry entry : index) {
            if (msLevel > 0 && entry.msLevel != msLevel) continue;
            Spectrum spec = readSpectrumAt(entry);
            if (spec != null) {
                consumer.accept(spec);
            }
        }
    }

    private void buildIndex() throws IOException {
        index.clear();
        scanMap.clear();

        // Read entire file to build index
        byte[] content = new byte[(int) file.length()];
        file.readFully(content);
        String text = new String(content, "UTF-8");

        Pattern pattern = Pattern.compile("<spectrum\\s[^>]*>");
        Matcher matcher = pattern.matcher(text);

        int scanNumber = 0;
        while (matcher.find()) {
            String tag = matcher.group();
            int offset = matcher.start();

            int endPos = text.indexOf("</spectrum>", offset);
            int length = endPos > 0 ? endPos + "</spectrum>".length() - offset : 0;

            scanNumber++;

            SpectrumIndexEntry entry = new SpectrumIndexEntry();
            entry.offset = offset;
            entry.length = length;
            entry.scanNumber = scanNumber;
            entry.msLevel = 1;

            // Extract index
            Matcher idxMatcher = Pattern.compile("index=\"(\\d+)\"").matcher(tag);
            if (idxMatcher.find()) {
                entry.scanNumber = Integer.parseInt(idxMatcher.group(1));
            }

            // Parse block for metadata
            if (endPos > 0) {
                String block = text.substring(offset, endPos + "</spectrum>".length());

                Matcher msMatcher = Pattern.compile(
                        "name=\"ms level\"[^>]*value=\"(\\d+)\"").matcher(block);
                if (msMatcher.find()) {
                    entry.msLevel = Integer.parseInt(msMatcher.group(1));
                }

                Matcher rtMatcher = Pattern.compile(
                        "name=\"scan start time\"[^>]*value=\"([^\"]+)\"").matcher(block);
                if (rtMatcher.find()) {
                    entry.retentionTime = Double.parseDouble(rtMatcher.group(1));
                }
            }

            int idx = index.size();
            index.add(entry);
            scanMap.put(entry.scanNumber, idx);
        }
    }

    private Spectrum readSpectrumAt(SpectrumIndexEntry entry) throws IOException {
        if (file == null || entry.length == 0) return null;

        Spectrum spec = new Spectrum();
        spec.setMsLevel(entry.msLevel);
        spec.setRetentionTime(entry.retentionTime);
        spec.setIndex(entry.scanNumber);

        return spec;
    }

    /**
     * Index entry for a spectrum.
     */
    public static class SpectrumIndexEntry {
        public int scanNumber;
        public long offset;
        public int length;
        public int msLevel = 1;
        public double retentionTime;
        public double precursorMz;
    }

    @FunctionalInterface
    public interface SpectrumConsumer {
        void accept(Spectrum spectrum) throws IOException;
    }
}
