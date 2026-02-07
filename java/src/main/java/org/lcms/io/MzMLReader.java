package org.lcms.io;

import org.lcms.core.*;

import javax.xml.parsers.*;
import org.xml.sax.*;
import org.xml.sax.helpers.*;
import java.io.*;
import java.util.*;
import java.util.zip.*;

/**
 * Reader for mzML format files.
 */
public class MzMLReader {

    public MSExperiment read(String filename) throws Exception {
        SAXParserFactory factory = SAXParserFactory.newInstance();
        factory.setNamespaceAware(true);
        SAXParser parser = factory.newSAXParser();

        MzMLHandler handler = new MzMLHandler();
        parser.parse(new File(filename), handler);

        MSExperiment exp = handler.getExperiment();
        exp.setSourceFile(filename);
        return exp;
    }

    private static class MzMLHandler extends DefaultHandler {
        private MSExperiment experiment;
        private Spectrum currentSpectrum;
        private Chromatogram currentChromatogram;
        private StringBuilder characterBuffer;
        private boolean inSpectrum;
        private boolean inChromatogram;
        private boolean inBinaryDataArray;

        // Binary data parsing
        private boolean isMzArray;
        private boolean isIntensityArray;
        private boolean isTimeArray;
        private boolean is64Bit;
        private boolean isCompressed;

        private List<double[]> binaryArrays;

        public MzMLHandler() {
            this.experiment = new MSExperiment();
            this.characterBuffer = new StringBuilder();
            this.binaryArrays = new ArrayList<>();
        }

        public MSExperiment getExperiment() {
            return experiment;
        }

        @Override
        public void startElement(String uri, String localName, String qName,
                                Attributes attributes) throws SAXException {
            characterBuffer.setLength(0);

            switch (localName) {
                case "spectrum":
                    inSpectrum = true;
                    currentSpectrum = new Spectrum();
                    currentSpectrum.setNativeId(attributes.getValue("id"));
                    binaryArrays.clear();
                    break;

                case "chromatogram":
                    inChromatogram = true;
                    currentChromatogram = new Chromatogram();
                    currentChromatogram.setNativeId(attributes.getValue("id"));
                    binaryArrays.clear();
                    break;

                case "binaryDataArray":
                    inBinaryDataArray = true;
                    isMzArray = false;
                    isIntensityArray = false;
                    isTimeArray = false;
                    is64Bit = true;
                    isCompressed = false;
                    break;

                case "cvParam":
                    handleCvParam(attributes);
                    break;

                case "precursor":
                    if (currentSpectrum != null) {
                        currentSpectrum.addPrecursor(new Precursor());
                    }
                    break;
            }
        }

        @Override
        public void endElement(String uri, String localName, String qName)
                throws SAXException {
            switch (localName) {
                case "spectrum":
                    if (binaryArrays.size() >= 2) {
                        currentSpectrum.setData(binaryArrays.get(0), binaryArrays.get(1));
                    }
                    experiment.addSpectrum(currentSpectrum);
                    currentSpectrum = null;
                    inSpectrum = false;
                    break;

                case "chromatogram":
                    if (binaryArrays.size() >= 2) {
                        currentChromatogram.setData(binaryArrays.get(0), binaryArrays.get(1));
                    }
                    experiment.addChromatogram(currentChromatogram);
                    currentChromatogram = null;
                    inChromatogram = false;
                    break;

                case "binaryDataArray":
                    inBinaryDataArray = false;
                    break;

                case "binary":
                    if (inBinaryDataArray) {
                        String base64Data = characterBuffer.toString().replaceAll("\\s", "");
                        if (!base64Data.isEmpty()) {
                            try {
                                double[] data = decodeBinary(base64Data, is64Bit, isCompressed);
                                binaryArrays.add(data);
                            } catch (Exception e) {
                                // Skip invalid binary data
                            }
                        }
                    }
                    break;
            }
        }

        @Override
        public void characters(char[] ch, int start, int length) {
            characterBuffer.append(ch, start, length);
        }

        private void handleCvParam(Attributes attrs) {
            String accession = attrs.getValue("accession");
            String value = attrs.getValue("value");
            String unitAccession = attrs.getValue("unitAccession");

            if (accession == null) return;

            if (inBinaryDataArray) {
                switch (accession) {
                    case "MS:1000514": isMzArray = true; break;
                    case "MS:1000515": isIntensityArray = true; break;
                    case "MS:1000595": isTimeArray = true; break;
                    case "MS:1000523": is64Bit = true; break;
                    case "MS:1000521": is64Bit = false; break;
                    case "MS:1000574": isCompressed = true; break;
                }
            } else if (currentSpectrum != null) {
                switch (accession) {
                    case "MS:1000511": // MS level
                        currentSpectrum.setMsLevel(Integer.parseInt(value));
                        break;
                    case "MS:1000128": // Profile
                        currentSpectrum.setType(SpectrumType.PROFILE);
                        break;
                    case "MS:1000127": // Centroid
                        currentSpectrum.setType(SpectrumType.CENTROID);
                        break;
                    case "MS:1000130": // Positive
                        currentSpectrum.setPolarity(Polarity.POSITIVE);
                        break;
                    case "MS:1000129": // Negative
                        currentSpectrum.setPolarity(Polarity.NEGATIVE);
                        break;
                    case "MS:1000016": // Scan start time
                        double rt = Double.parseDouble(value);
                        if ("UO:0000031".equals(unitAccession)) {
                            rt *= 60; // Convert minutes to seconds
                        }
                        currentSpectrum.setRetentionTime(rt);
                        break;
                    case "MS:1000744": // Selected ion m/z
                        if (!currentSpectrum.getPrecursors().isEmpty()) {
                            Precursor p = currentSpectrum.getPrecursors()
                                .get(currentSpectrum.getPrecursors().size() - 1);
                            p.setMz(Double.parseDouble(value));
                        }
                        break;
                    case "MS:1000041": // Charge state
                        if (!currentSpectrum.getPrecursors().isEmpty()) {
                            Precursor p = currentSpectrum.getPrecursors()
                                .get(currentSpectrum.getPrecursors().size() - 1);
                            p.setCharge(Integer.parseInt(value));
                        }
                        break;
                }
            } else if (currentChromatogram != null) {
                switch (accession) {
                    case "MS:1000235": // TIC
                        currentChromatogram.setType(ChromatogramType.TIC);
                        break;
                    case "MS:1000628": // BPC
                        currentChromatogram.setType(ChromatogramType.BPC);
                        break;
                    case "MS:1001473": // SRM
                        currentChromatogram.setType(ChromatogramType.SRM);
                        break;
                }
            }
        }

        private double[] decodeBinary(String base64, boolean is64Bit, boolean isCompressed)
                throws Exception {
            byte[] decoded = Base64.getDecoder().decode(base64);

            if (isCompressed) {
                decoded = decompress(decoded);
            }

            int elementSize = is64Bit ? 8 : 4;
            int count = decoded.length / elementSize;
            double[] result = new double[count];

            for (int i = 0; i < count; i++) {
                if (is64Bit) {
                    long bits = 0;
                    for (int j = 0; j < 8; j++) {
                        bits |= ((long) (decoded[i * 8 + j] & 0xFF)) << (j * 8);
                    }
                    result[i] = Double.longBitsToDouble(bits);
                } else {
                    int bits = 0;
                    for (int j = 0; j < 4; j++) {
                        bits |= (decoded[i * 4 + j] & 0xFF) << (j * 8);
                    }
                    result[i] = Float.intBitsToFloat(bits);
                }
            }

            return result;
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
