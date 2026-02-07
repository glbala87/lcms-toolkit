package org.lcms.io;

import org.lcms.core.*;

import javax.xml.parsers.*;
import org.xml.sax.*;
import org.xml.sax.helpers.*;
import java.io.*;
import java.util.*;
import java.util.zip.*;

/**
 * Reader for mzXML format files.
 */
public class MzXMLReader {

    public MSExperiment read(String filename) throws Exception {
        SAXParserFactory factory = SAXParserFactory.newInstance();
        factory.setNamespaceAware(true);
        SAXParser parser = factory.newSAXParser();

        MzXMLHandler handler = new MzXMLHandler();
        parser.parse(new File(filename), handler);

        MSExperiment exp = handler.getExperiment();
        exp.setSourceFile(filename);
        return exp;
    }

    private static class MzXMLHandler extends DefaultHandler {
        private MSExperiment experiment;
        private Spectrum currentSpectrum;
        private StringBuilder characterBuffer;
        private boolean inScan;
        private boolean inPeaks;

        // Peak data parsing
        private int precision;
        private String byteOrder;
        private String compressionType;
        private String pairOrder;

        public MzXMLHandler() {
            this.experiment = new MSExperiment();
            this.characterBuffer = new StringBuilder();
        }

        public MSExperiment getExperiment() {
            return experiment;
        }

        @Override
        public void startElement(String uri, String localName, String qName,
                                Attributes attributes) throws SAXException {
            characterBuffer.setLength(0);

            switch (localName) {
                case "scan":
                    inScan = true;
                    currentSpectrum = new Spectrum();
                    parseScanAttributes(attributes);
                    break;

                case "peaks":
                    inPeaks = true;
                    precision = Integer.parseInt(attributes.getValue("precision") != null ?
                        attributes.getValue("precision") : "32");
                    byteOrder = attributes.getValue("byteOrder");
                    compressionType = attributes.getValue("compressionType");
                    pairOrder = attributes.getValue("pairOrder");
                    break;

                case "precursorMz":
                    if (currentSpectrum != null) {
                        Precursor precursor = new Precursor();

                        String precursorIntensity = attributes.getValue("precursorIntensity");
                        if (precursorIntensity != null) {
                            precursor.setIntensity(Double.parseDouble(precursorIntensity));
                        }

                        String precursorCharge = attributes.getValue("precursorCharge");
                        if (precursorCharge != null) {
                            precursor.setCharge(Integer.parseInt(precursorCharge));
                        }

                        String activationMethod = attributes.getValue("activationMethod");
                        if (activationMethod != null) {
                            try {
                                precursor.setActivationMethod(
                                    ActivationMethod.valueOf(activationMethod));
                            } catch (IllegalArgumentException e) {
                                precursor.setActivationMethod(ActivationMethod.UNKNOWN);
                            }
                        }

                        currentSpectrum.addPrecursor(precursor);
                    }
                    break;

                case "msInstrument":
                    break;

                case "msModel":
                    if (attributes.getValue("value") != null) {
                        experiment.setInstrumentModel(attributes.getValue("value"));
                    }
                    break;
            }
        }

        @Override
        public void endElement(String uri, String localName, String qName)
                throws SAXException {
            switch (localName) {
                case "scan":
                    if (currentSpectrum != null) {
                        experiment.addSpectrum(currentSpectrum);
                    }
                    currentSpectrum = null;
                    inScan = false;
                    break;

                case "peaks":
                    if (currentSpectrum != null && inPeaks) {
                        String base64Data = characterBuffer.toString().replaceAll("\\s", "");
                        if (!base64Data.isEmpty()) {
                            try {
                                parseAndSetPeaks(base64Data);
                            } catch (Exception e) {
                                // Skip invalid peak data
                            }
                        }
                    }
                    inPeaks = false;
                    break;

                case "precursorMz":
                    if (currentSpectrum != null && !currentSpectrum.getPrecursors().isEmpty()) {
                        String mzStr = characterBuffer.toString().trim();
                        if (!mzStr.isEmpty()) {
                            Precursor p = currentSpectrum.getPrecursors()
                                .get(currentSpectrum.getPrecursors().size() - 1);
                            p.setMz(Double.parseDouble(mzStr));
                        }
                    }
                    break;
            }
        }

        @Override
        public void characters(char[] ch, int start, int length) {
            characterBuffer.append(ch, start, length);
        }

        private void parseScanAttributes(Attributes attrs) {
            String msLevel = attrs.getValue("msLevel");
            if (msLevel != null) {
                currentSpectrum.setMsLevel(Integer.parseInt(msLevel));
            }

            String retentionTime = attrs.getValue("retentionTime");
            if (retentionTime != null) {
                double rt = parseRetentionTime(retentionTime);
                currentSpectrum.setRetentionTime(rt);
            }

            String polarity = attrs.getValue("polarity");
            if (polarity != null) {
                if (polarity.equals("+") || polarity.equalsIgnoreCase("positive")) {
                    currentSpectrum.setPolarity(Polarity.POSITIVE);
                } else if (polarity.equals("-") || polarity.equalsIgnoreCase("negative")) {
                    currentSpectrum.setPolarity(Polarity.NEGATIVE);
                }
            }

            String centroided = attrs.getValue("centroided");
            if ("1".equals(centroided) || "true".equalsIgnoreCase(centroided)) {
                currentSpectrum.setType(SpectrumType.CENTROID);
            } else {
                currentSpectrum.setType(SpectrumType.PROFILE);
            }
        }

        private double parseRetentionTime(String rtStr) {
            // Format: PT123.456S or PT1.5M or just a number
            if (rtStr.startsWith("PT") && rtStr.endsWith("S")) {
                return Double.parseDouble(rtStr.substring(2, rtStr.length() - 1));
            } else if (rtStr.startsWith("PT") && rtStr.endsWith("M")) {
                return Double.parseDouble(rtStr.substring(2, rtStr.length() - 1)) * 60;
            } else {
                return Double.parseDouble(rtStr);
            }
        }

        private void parseAndSetPeaks(String base64Data) throws Exception {
            byte[] decoded = Base64.getDecoder().decode(base64Data);

            // Handle compression
            if ("zlib".equals(compressionType)) {
                decoded = decompress(decoded);
            }

            // Parse floats/doubles
            boolean is64Bit = (precision == 64);
            boolean isLittleEndian = !"network".equals(byteOrder);
            boolean mzFirst = !"int-mz".equals(pairOrder);

            int elementSize = is64Bit ? 8 : 4;
            int pairSize = elementSize * 2;
            int pairCount = decoded.length / pairSize;

            double[] mz = new double[pairCount];
            double[] intensity = new double[pairCount];

            for (int i = 0; i < pairCount; i++) {
                double val1, val2;

                if (is64Bit) {
                    val1 = readDouble(decoded, i * pairSize, isLittleEndian);
                    val2 = readDouble(decoded, i * pairSize + 8, isLittleEndian);
                } else {
                    val1 = readFloat(decoded, i * pairSize, isLittleEndian);
                    val2 = readFloat(decoded, i * pairSize + 4, isLittleEndian);
                }

                if (mzFirst) {
                    mz[i] = val1;
                    intensity[i] = val2;
                } else {
                    intensity[i] = val1;
                    mz[i] = val2;
                }
            }

            currentSpectrum.setData(mz, intensity);
        }

        private double readDouble(byte[] data, int offset, boolean littleEndian) {
            long bits = 0;
            if (littleEndian) {
                for (int i = 0; i < 8; i++) {
                    bits |= ((long) (data[offset + i] & 0xFF)) << (i * 8);
                }
            } else {
                for (int i = 0; i < 8; i++) {
                    bits |= ((long) (data[offset + i] & 0xFF)) << ((7 - i) * 8);
                }
            }
            return Double.longBitsToDouble(bits);
        }

        private float readFloat(byte[] data, int offset, boolean littleEndian) {
            int bits = 0;
            if (littleEndian) {
                for (int i = 0; i < 4; i++) {
                    bits |= (data[offset + i] & 0xFF) << (i * 8);
                }
            } else {
                for (int i = 0; i < 4; i++) {
                    bits |= (data[offset + i] & 0xFF) << ((3 - i) * 8);
                }
            }
            return Float.intBitsToFloat(bits);
        }

        private byte[] decompress(byte[] data) throws Exception {
            Inflater inflater = new Inflater();
            inflater.setInput(data);

            ByteArrayOutputStream outputStream = new ByteArrayOutputStream(data.length * 4);
            byte[] buffer = new byte[1024];

            while (!inflater.finished()) {
                int count = inflater.inflate(buffer);
                outputStream.write(buffer, 0, count);
            }

            outputStream.close();
            inflater.end();
            return outputStream.toByteArray();
        }
    }
}
