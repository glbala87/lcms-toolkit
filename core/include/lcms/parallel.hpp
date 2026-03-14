#pragma once

#include "spectrum.hpp"
#include "peak.hpp"
#include "ms_experiment.hpp"
#include <vector>
#include <functional>
#include <thread>
#include <mutex>
#include <future>
#include <queue>

namespace lcms {

/**
 * @brief Simple thread pool for parallel processing.
 */
class ThreadPool {
public:
    explicit ThreadPool(size_t num_threads = 0)
        : stop_(false) {
        if (num_threads == 0) {
            num_threads = std::thread::hardware_concurrency();
            if (num_threads == 0) num_threads = 4;
        }

        for (size_t i = 0; i < num_threads; ++i) {
            workers_.emplace_back([this] {
                while (true) {
                    std::function<void()> task;
                    {
                        std::unique_lock<std::mutex> lock(mutex_);
                        cv_.wait(lock, [this] {
                            return stop_ || !tasks_.empty();
                        });
                        if (stop_ && tasks_.empty()) return;
                        task = std::move(tasks_.front());
                        tasks_.pop();
                    }
                    task();
                }
            });
        }
    }

    ~ThreadPool() {
        {
            std::unique_lock<std::mutex> lock(mutex_);
            stop_ = true;
        }
        cv_.notify_all();
        for (auto& worker : workers_) {
            if (worker.joinable()) worker.join();
        }
    }

    /// Submit a task and get a future
    template<typename F, typename... Args>
    auto submit(F&& f, Args&&... args)
        -> std::future<decltype(f(args...))> {
        using ReturnType = decltype(f(args...));

        auto task = std::make_shared<std::packaged_task<ReturnType()>>(
            std::bind(std::forward<F>(f), std::forward<Args>(args)...)
        );

        std::future<ReturnType> result = task->get_future();
        {
            std::unique_lock<std::mutex> lock(mutex_);
            tasks_.emplace([task]() { (*task)(); });
        }
        cv_.notify_one();
        return result;
    }

    [[nodiscard]] size_t numThreads() const { return workers_.size(); }

private:
    std::vector<std::thread> workers_;
    std::queue<std::function<void()>> tasks_;
    std::mutex mutex_;
    std::condition_variable cv_;
    bool stop_;
};

/**
 * @brief Process spectra in parallel.
 *
 * @param spectra Input spectra
 * @param func Function to apply to each spectrum
 * @param num_threads Number of threads (0 = auto)
 * @return Results from each spectrum
 */
template<typename ResultType>
std::vector<ResultType> parallelProcessSpectra(
    const std::vector<Spectrum>& spectra,
    std::function<ResultType(const Spectrum&)> func,
    size_t num_threads = 0) {

    ThreadPool pool(num_threads);
    std::vector<std::future<ResultType>> futures;
    futures.reserve(spectra.size());

    for (const auto& spectrum : spectra) {
        futures.push_back(pool.submit(func, std::cref(spectrum)));
    }

    std::vector<ResultType> results;
    results.reserve(futures.size());
    for (auto& future : futures) {
        results.push_back(future.get());
    }

    return results;
}

/**
 * @brief Parallel peak picking across multiple spectra.
 */
PeakList parallelPeakPicking(
    const std::vector<Spectrum>& spectra,
    double min_snr = 3.0,
    size_t num_threads = 0);

} // namespace lcms
