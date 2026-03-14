#include "lcms/io/indexed_reader.hpp"
#include <sstream>
#include <regex>
#include <algorithm>
#include <set>

namespace lcms {
namespace io {

// =========================================================================
// FileIndex
// =========================================================================

void FileIndex::addEntry(const SpectrumIndexEntry& entry) {
    size_t idx = entries_.size();
    entries_.push_back(entry);
    scan_map_[entry.scan_number] = idx;
}

const SpectrumIndexEntry* FileIndex::getByIndex(size_t index) const {
    if (index < entries_.size()) return &entries_[index];
    return nullptr;
}

const SpectrumIndexEntry* FileIndex::getByScan(Index scan) const {
    auto it = scan_map_.find(scan);
    if (it != scan_map_.end()) return &entries_[it->second];
    return nullptr;
}

const SpectrumIndexEntry* FileIndex::getByRT(RetentionTime rt, double tolerance) const {
    const SpectrumIndexEntry* best = nullptr;
    double best_diff = tolerance + 1.0;

    for (const auto& entry : entries_) {
        double diff = std::abs(entry.rt - rt);
        if (diff <= tolerance && diff < best_diff) {
            best = &entry;
            best_diff = diff;
        }
    }

    return best;
}

RTRange FileIndex::rtRange() const {
    RTRange range;
    for (const auto& entry : entries_) {
        range.extend(entry.rt);
    }
    return range;
}

std::vector<MSLevel> FileIndex::msLevels() const {
    std::set<MSLevel> levels;
    for (const auto& entry : entries_) {
        levels.insert(entry.ms_level);
    }
    return std::vector<MSLevel>(levels.begin(), levels.end());
}

void FileIndex::save(const std::string& filepath) const {
    std::ofstream file(filepath);
    file << entries_.size() << "\n";
    for (const auto& entry : entries_) {
        file << entry.scan_number << " "
             << entry.offset << " "
             << entry.length << " "
             << static_cast<int>(entry.ms_level) << " "
             << entry.rt << " "
             << entry.precursor_mz << "\n";
    }
}

FileIndex FileIndex::load(const std::string& filepath) {
    FileIndex index;
    std::ifstream file(filepath);
    if (!file.is_open()) return index;

    size_t count;
    file >> count;

    for (size_t i = 0; i < count; ++i) {
        SpectrumIndexEntry entry;
        int ms_level;
        file >> entry.scan_number >> entry.offset >> entry.length
             >> ms_level >> entry.rt >> entry.precursor_mz;
        entry.ms_level = static_cast<MSLevel>(ms_level);
        index.addEntry(entry);
    }

    return index;
}

FileIndex FileIndex::buildFromMzML(const std::string& filepath) {
    FileIndex index;

    std::ifstream file(filepath);
    if (!file.is_open()) return index;

    std::string content((std::istreambuf_iterator<char>(file)),
                         std::istreambuf_iterator<char>());

    std::regex spectrum_re("<spectrum\\s[^>]*>");
    Index scan_number = 0;

    auto begin = std::sregex_iterator(content.begin(), content.end(), spectrum_re);
    auto end = std::sregex_iterator();

    for (auto it = begin; it != end; ++it) {
        std::string tag = it->str();
        auto offset = static_cast<std::streamoff>(it->position());

        // Find end of spectrum
        std::string end_tag = "</spectrum>";
        auto end_pos = content.find(end_tag, offset);
        std::streamoff length = (end_pos != std::string::npos)
            ? static_cast<std::streamoff>(end_pos + end_tag.size() - offset) : 0;

        scan_number++;

        SpectrumIndexEntry entry;
        entry.scan_number = scan_number;
        entry.offset = offset;
        entry.length = length;
        entry.ms_level = 1;

        // Extract index attribute
        std::regex idx_re("index=\"(\\d+)\"");
        std::smatch idx_match;
        if (std::regex_search(tag, idx_match, idx_re)) {
            entry.scan_number = std::stoull(idx_match[1].str());
        }

        // Parse block for metadata
        if (end_pos != std::string::npos) {
            std::string block = content.substr(offset, end_pos + end_tag.size() - offset);

            std::regex ms_re("name=\"ms level\"[^>]*value=\"(\\d+)\"");
            std::smatch ms_match;
            if (std::regex_search(block, ms_match, ms_re)) {
                entry.ms_level = static_cast<MSLevel>(std::stoi(ms_match[1].str()));
            }

            std::regex rt_re("name=\"scan start time\"[^>]*value=\"([^\"]+)\"");
            std::smatch rt_match;
            if (std::regex_search(block, rt_match, rt_re)) {
                entry.rt = std::stod(rt_match[1].str());
            }

            std::regex pre_re("name=\"selected ion m/z\"[^>]*value=\"([^\"]+)\"");
            std::smatch pre_match;
            if (std::regex_search(block, pre_match, pre_re)) {
                entry.precursor_mz = std::stod(pre_match[1].str());
            }
        }

        index.addEntry(entry);
    }

    return index;
}

FileIndex FileIndex::buildFromMzXML(const std::string& filepath) {
    FileIndex index;

    std::ifstream file(filepath);
    if (!file.is_open()) return index;

    std::string content((std::istreambuf_iterator<char>(file)),
                         std::istreambuf_iterator<char>());

    std::regex scan_re("<scan\\s[^>]*>");
    auto begin = std::sregex_iterator(content.begin(), content.end(), scan_re);
    auto end = std::sregex_iterator();

    for (auto it = begin; it != end; ++it) {
        std::string tag = it->str();
        auto offset = static_cast<std::streamoff>(it->position());

        std::string end_tag = "</scan>";
        auto end_pos = content.find(end_tag, offset);
        std::streamoff length = (end_pos != std::string::npos)
            ? static_cast<std::streamoff>(end_pos + end_tag.size() - offset) : 0;

        SpectrumIndexEntry entry;
        entry.offset = offset;
        entry.length = length;

        std::regex num_re("num=\"(\\d+)\"");
        std::smatch num_match;
        if (std::regex_search(tag, num_match, num_re)) {
            entry.scan_number = std::stoull(num_match[1].str());
        }

        std::regex level_re("msLevel=\"(\\d+)\"");
        std::smatch level_match;
        if (std::regex_search(tag, level_match, level_re)) {
            entry.ms_level = static_cast<MSLevel>(std::stoi(level_match[1].str()));
        }

        std::regex rt_re("retentionTime=\"PT([\\d.]+)S\"");
        std::smatch rt_match;
        if (std::regex_search(tag, rt_match, rt_re)) {
            entry.rt = std::stod(rt_match[1].str());
        }

        index.addEntry(entry);
    }

    return index;
}

// =========================================================================
// IndexedMzMLReader
// =========================================================================

IndexedMzMLReader::IndexedMzMLReader(const std::string& filepath)
    : filepath_(filepath) {}

IndexedMzMLReader::~IndexedMzMLReader() {
    close();
}

bool IndexedMzMLReader::open() {
    file_ = std::make_unique<std::ifstream>(filepath_, std::ios::binary);
    if (!file_->is_open()) return false;
    is_open_ = true;

    // Try loading cached index
    std::string idx_path = filepath_ + ".idx";
    std::ifstream idx_file(idx_path);
    if (idx_file.good()) {
        idx_file.close();
        index_ = FileIndex::load(idx_path);
    } else {
        index_ = FileIndex::buildFromMzML(filepath_);
    }

    return true;
}

void IndexedMzMLReader::close() {
    if (file_) {
        file_->close();
        file_.reset();
    }
    is_open_ = false;
}

size_t IndexedMzMLReader::spectrumCount() const {
    return index_.size();
}

std::optional<Spectrum> IndexedMzMLReader::getSpectrum(size_t idx) {
    auto entry = index_.getByIndex(idx);
    if (!entry) return std::nullopt;
    return readSpectrumAt(*entry);
}

std::optional<Spectrum> IndexedMzMLReader::getSpectrumByScan(Index scan) {
    auto entry = index_.getByScan(scan);
    if (!entry) return std::nullopt;
    return readSpectrumAt(*entry);
}

std::optional<Spectrum> IndexedMzMLReader::getSpectrumByRT(
    RetentionTime rt, double tolerance) {
    auto entry = index_.getByRT(rt, tolerance);
    if (!entry) return std::nullopt;
    return readSpectrumAt(*entry);
}

void IndexedMzMLReader::saveIndex(const std::string& filepath) const {
    std::string path = filepath.empty() ? (filepath_ + ".idx") : filepath;
    index_.save(path);
}

std::optional<Spectrum> IndexedMzMLReader::readSpectrumAt(
    const SpectrumIndexEntry& entry) {
    if (!is_open_ || !file_ || entry.length == 0) return std::nullopt;

    file_->seekg(entry.offset);
    std::string data(entry.length, '\0');
    file_->read(&data[0], entry.length);

    Spectrum spectrum;
    spectrum.setMsLevel(entry.ms_level);
    spectrum.setRt(entry.rt);
    spectrum.setIndex(entry.scan_number);

    return spectrum;
}

// =========================================================================
// IndexedMzXMLReader
// =========================================================================

IndexedMzXMLReader::IndexedMzXMLReader(const std::string& filepath)
    : filepath_(filepath) {}

IndexedMzXMLReader::~IndexedMzXMLReader() {
    close();
}

bool IndexedMzXMLReader::open() {
    file_ = std::make_unique<std::ifstream>(filepath_, std::ios::binary);
    if (!file_->is_open()) return false;
    is_open_ = true;

    std::string idx_path = filepath_ + ".idx";
    std::ifstream idx_file(idx_path);
    if (idx_file.good()) {
        idx_file.close();
        index_ = FileIndex::load(idx_path);
    } else {
        index_ = FileIndex::buildFromMzXML(filepath_);
    }

    return true;
}

void IndexedMzXMLReader::close() {
    if (file_) {
        file_->close();
        file_.reset();
    }
    is_open_ = false;
}

size_t IndexedMzXMLReader::spectrumCount() const {
    return index_.size();
}

std::optional<Spectrum> IndexedMzXMLReader::getSpectrum(size_t idx) {
    auto entry = index_.getByIndex(idx);
    if (!entry) return std::nullopt;
    return readSpectrumAt(*entry);
}

std::optional<Spectrum> IndexedMzXMLReader::getSpectrumByScan(Index scan) {
    auto entry = index_.getByScan(scan);
    if (!entry) return std::nullopt;
    return readSpectrumAt(*entry);
}

std::optional<Spectrum> IndexedMzXMLReader::getSpectrumByRT(
    RetentionTime rt, double tolerance) {
    auto entry = index_.getByRT(rt, tolerance);
    if (!entry) return std::nullopt;
    return readSpectrumAt(*entry);
}

void IndexedMzXMLReader::saveIndex(const std::string& filepath) const {
    std::string path = filepath.empty() ? (filepath_ + ".idx") : filepath;
    index_.save(path);
}

std::optional<Spectrum> IndexedMzXMLReader::readSpectrumAt(
    const SpectrumIndexEntry& entry) {
    if (!is_open_ || !file_ || entry.length == 0) return std::nullopt;

    file_->seekg(entry.offset);
    std::string data(entry.length, '\0');
    file_->read(&data[0], entry.length);

    Spectrum spectrum;
    spectrum.setMsLevel(entry.ms_level);
    spectrum.setRt(entry.rt);
    spectrum.setIndex(entry.scan_number);

    return spectrum;
}

} // namespace io
} // namespace lcms
