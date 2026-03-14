package org.lcms.tools;

import org.lcms.core.*;

import java.util.*;
import java.util.concurrent.*;
import java.util.function.Function;

/**
 * Parallel processor for spectra within a single experiment.
 */
public class ParallelProcessor {

    private final int numThreads;

    public ParallelProcessor() {
        this(Runtime.getRuntime().availableProcessors());
    }

    public ParallelProcessor(int numThreads) {
        this.numThreads = numThreads;
    }

    /**
     * Process spectra in parallel.
     *
     * @param spectra   List of spectra to process
     * @param processor Function to apply to each spectrum
     * @return List of results
     */
    public <T> List<T> processSpectra(List<Spectrum> spectra,
                                       Function<Spectrum, T> processor) {
        ExecutorService executor = Executors.newFixedThreadPool(numThreads);

        try {
            List<Future<T>> futures = new ArrayList<>();
            for (Spectrum spec : spectra) {
                futures.add(executor.submit(() -> processor.apply(spec)));
            }

            List<T> results = new ArrayList<>();
            for (Future<T> future : futures) {
                try {
                    results.add(future.get());
                } catch (Exception e) {
                    results.add(null);
                }
            }

            return results;
        } finally {
            executor.shutdown();
        }
    }

    /**
     * Process spectra in chunks for memory efficiency.
     *
     * @param spectra   List of spectra
     * @param chunkSize Number of spectra per chunk
     * @param processor Function to process a chunk of spectra
     * @return Aggregated results
     */
    public <T> List<T> processChunks(List<Spectrum> spectra, int chunkSize,
                                      Function<List<Spectrum>, T> processor) {
        List<T> results = new ArrayList<>();

        for (int i = 0; i < spectra.size(); i += chunkSize) {
            int end = Math.min(i + chunkSize, spectra.size());
            List<Spectrum> chunk = spectra.subList(i, end);
            results.add(processor.apply(chunk));
        }

        return results;
    }

    public int getNumThreads() { return numThreads; }
}
