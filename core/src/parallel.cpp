#include "lcms/parallel.hpp"
#include "lcms/algorithms/peak_picking.hpp"

namespace lcms {

PeakList parallelPeakPicking(
    const std::vector<Spectrum>& spectra,
    double min_snr,
    size_t num_threads) {

    algorithms::PeakPickingOptions options;
    options.min_snr = min_snr;

    auto pick_func = [&options](const Spectrum& spectrum) -> PeakList {
        algorithms::PeakPicker picker(options);
        return picker.pick(spectrum);
    };

    auto results = parallelProcessSpectra<PeakList>(spectra, pick_func, num_threads);

    // Merge all peak lists
    PeakList merged;
    for (auto& pl : results) {
        for (size_t i = 0; i < pl.size(); ++i) {
            merged.add(pl[i]);
        }
    }

    return merged;
}

} // namespace lcms
