package org.lcms.io;

import org.lcms.core.*;

import java.io.*;
import java.util.*;
import java.util.regex.*;

/**
 * Indexed mzXML reader for random access to spectra.
 */
public class IndexedMzXMLReader implements AutoCloseable {

    private final String filepath;
    private RandomAccessFile file;
    private List<SpectrumIndexEntry> index;
    private Map<Integer, Integer> scanMap;

    public IndexedMzXMLReader(String filepath) {
        this.filepath = filepath;
        this.index = new ArrayList<>();
        this.scanMap = new HashMap<>();
    }

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

    public Spectrum getSpectrum(int idx) throws IOException {
        if (idx < 0 || idx >= index.size()) return null;
        return readSpectrumAt(index.get(idx));
    }

    public Spectrum getSpectrumByScan(int scanNumber) throws IOException {
        Integer idx = scanMap.get(scanNumber);
        if (idx == null) return null;
        return readSpectrumAt(index.get(idx));
    }

    public void iterateSpectra(IndexedMzMLReader.SpectrumConsumer consumer) throws IOException {
        iterateSpectra(consumer, -1);
    }

    public void iterateSpectra(IndexedMzMLReader.SpectrumConsumer consumer,
                                int msLevel) throws IOException {
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

        byte[] content = new byte[(int) file.length()];
        file.readFully(content);
        String text = new String(content, "UTF-8");

        Pattern pattern = Pattern.compile("<scan\\s[^>]*>");
        Matcher matcher = pattern.matcher(text);

        while (matcher.find()) {
            String tag = matcher.group();
            int offset = matcher.start();

            int endPos = text.indexOf("</scan>", offset);
            int length = endPos > 0 ? endPos + "</scan>".length() - offset : 0;

            SpectrumIndexEntry entry = new SpectrumIndexEntry();
            entry.offset = offset;
            entry.length = length;

            Matcher numMatcher = Pattern.compile("num=\"(\\d+)\"").matcher(tag);
            if (numMatcher.find()) {
                entry.scanNumber = Integer.parseInt(numMatcher.group(1));
            }

            Matcher levelMatcher = Pattern.compile("msLevel=\"(\\d+)\"").matcher(tag);
            if (levelMatcher.find()) {
                entry.msLevel = Integer.parseInt(levelMatcher.group(1));
            }

            Matcher rtMatcher = Pattern.compile("retentionTime=\"PT([\\d.]+)S\"").matcher(tag);
            if (rtMatcher.find()) {
                entry.retentionTime = Double.parseDouble(rtMatcher.group(1));
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

    public static class SpectrumIndexEntry {
        public int scanNumber;
        public long offset;
        public int length;
        public int msLevel = 1;
        public double retentionTime;
        public double precursorMz;
    }
}
