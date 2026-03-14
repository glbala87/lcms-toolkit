package org.lcms.tools;

import org.lcms.core.*;

import java.nio.file.Path;
import java.util.*;
import java.util.concurrent.*;
import java.util.function.Function;

/**
 * Parallel batch processor for LC-MS files.
 */
public class BatchProcessor {

    private final int numWorkers;
    private final ExecutorService executor;

    public BatchProcessor() {
        this(Runtime.getRuntime().availableProcessors());
    }

    public BatchProcessor(int numWorkers) {
        this.numWorkers = numWorkers;
        this.executor = Executors.newFixedThreadPool(numWorkers);
    }

    /**
     * Process multiple files in parallel.
     *
     * @param files     List of file paths
     * @param processor Function that processes a path and returns a result
     * @return Map of file path to result
     */
    public <T> Map<String, T> processFiles(List<Path> files,
                                             Function<Path, T> processor) {
        Map<String, Future<T>> futures = new LinkedHashMap<>();
        for (Path file : files) {
            futures.put(file.toString(), executor.submit(() -> processor.apply(file)));
        }

        Map<String, T> results = new LinkedHashMap<>();
        for (Map.Entry<String, Future<T>> entry : futures.entrySet()) {
            try {
                results.put(entry.getKey(), entry.getValue().get());
            } catch (Exception e) {
                System.err.println("Error processing " + entry.getKey() + ": " + e.getMessage());
            }
        }

        return results;
    }

    /**
     * Apply a function to items in parallel.
     */
    public <T, R> List<R> map(List<T> items, Function<T, R> func) {
        List<Future<R>> futures = new ArrayList<>();
        for (T item : items) {
            futures.add(executor.submit(() -> func.apply(item)));
        }

        List<R> results = new ArrayList<>();
        for (Future<R> future : futures) {
            try {
                results.add(future.get());
            } catch (Exception e) {
                results.add(null);
            }
        }

        return results;
    }

    /**
     * Shutdown the thread pool.
     */
    public void shutdown() {
        executor.shutdown();
        try {
            if (!executor.awaitTermination(60, TimeUnit.SECONDS)) {
                executor.shutdownNow();
            }
        } catch (InterruptedException e) {
            executor.shutdownNow();
            Thread.currentThread().interrupt();
        }
    }

    public int getNumWorkers() { return numWorkers; }
}
