#include <future>
#include "CSVPreprocessor.h"

#include "Utils.h"
#include "string_table_io/StringTableWriter.h"
#include "fingerprint_table_io/FingerprintTableWriter.h"
#include "id_to_string_io/IdToStringWriter.h"
#include "properties_table_io/PropertiesTableWriter.h"
#include "PropertiesFilter.h"
#include "Profiling.h"

using namespace std;
using namespace qtr;

namespace {
    std::string strip(const std::string& input) {
        size_t start = 0;
        size_t end = input.length();
        auto is_special_char = [](char c) {
            return std::isspace(c) || c == '\n' || c == '\t' || c == '\r';
        };
        while (start < end && is_special_char(input[start])) {
            ++start;
        }
        while (end > start && is_special_char(input[end - 1])) {
            --end;
        }
        return input.substr(start, end - start);
    }

    std::vector<std::string> splitString(const std::string &str) {
        static const char delimiter = '\t';
        std::vector<std::string> substrings;
        std::size_t start = 0;
        std::size_t end = str.find(delimiter);
        while (end != std::string::npos) {
            substrings.push_back(strip(str.substr(start, end - start)));
            start = end + 1;
            end = str.find(delimiter, start);
        }
        substrings.push_back(strip(str.substr(start, end)));
        return substrings;
    }

    uint64_t
    parseCSV(const std::filesystem::path &csvFilePath, const PreprocessingArgs &args, std::atomic_uint64_t &counter,
             std::mutex &mutex) {
        std::filesystem::path smilesTablePath =
                args.smilesTables() / (csvFilePath.stem().string() + stringTableExtension);
        std::filesystem::path fingerprintTablePath =
                args.fingerprintTables() / (csvFilePath.stem().string() + fingerprintTableExtension);
        std::filesystem::path idToStringTablePath =
                args.idToStringTables() / (csvFilePath.stem().string() + ".csv");
        std::filesystem::path propertyTablesPath = args.propertyTables() / csvFilePath.stem();

        std::filesystem::create_directory(smilesTablePath.parent_path());
        std::filesystem::create_directory(fingerprintTablePath.parent_path());
        std::filesystem::create_directory(idToStringTablePath.parent_path());
        if (args.properties())
            std::filesystem::create_directory(propertyTablesPath.parent_path());

        StringTableWriter smilesTableWriter(smilesTablePath);
        FingerprintTableWriter fingerprintTableWriter(fingerprintTablePath);
        IdToStringWriter idToStringWriter(idToStringTablePath);
        unique_ptr<PropertiesTableWriter> propertiesTableWriter = nullptr;
        if (args.properties())
            propertiesTableWriter = make_unique<PropertiesTableWriter>(propertyTablesPath);

        std::ifstream in(csvFilePath);
        std::string line;
        uint64_t skipped = 0, processed = 0;
        uint64_t fingerprintLength = 0;
        while (std::getline(in, line)) {
            auto lineElements = splitString(line);
            size_t expectedElements = 2 + size_t(args.properties()) * PropertiesFilter::Properties().size() +
                                      size_t(args.fingerprintType() == qtr::FingerprintType::Custom);
            if (lineElements.size() != expectedElements) {
                LOG(WARNING) << "Skip line with invalid (" << lineElements.size()
                             << ") number of arguments elements (expected" << expectedElements << ": " << line;
                skipped++;
                continue;
            }

            std::string strId = lineElements[0];
            std::string smiles = lineElements[1];
            PropertiesFilter::Properties properties{};
            if (args.properties()) {
                try {
                    for (size_t i = 0; i < properties.size(); i++) {
                        properties[i] = std::stof(lineElements[i + 2]);
                    }
                }
                catch (const std::invalid_argument &e) {
                    LOG(WARNING) << "Skip line with invalid property: " << line;
                    skipped++;
                    continue;
                }
            }

            try {
                Fingerprint fingerprint = [&]() {
                    if (args.fingerprintType() == qtr::FingerprintType::Custom) {
                        size_t fingerprintIndex = 2 + size_t(args.properties()) * PropertiesFilter::Properties().size();
                        std::string fingerprintStr = lineElements[fingerprintIndex];
                        return Fingerprint(fingerprintStr);
                    } else if (args.fingerprintType() == qtr::FingerprintType::Indigo) {
                        return indigoFingerprintFromSmiles(smiles);
                    } else if (args.fingerprintType() == qtr::FingerprintType::RDKit) {
                        return rdkitFingerprintFromSmiles(smiles);
                    } else {
                        throw std::logic_error("Undefined fingerprint type specified");
                    }
                }();

                if (fingerprintLength == 0)
                    fingerprintLength = fingerprint.size();
                else
                    assert(fingerprintLength == fingerprint.size());

                {
                    std::lock_guard lock(mutex);
                    uint64_t id = counter++;
                    smilesTableWriter << std::make_pair(id, smiles);
                    fingerprintTableWriter << std::make_pair(id, fingerprint);
                    if (args.molIdType() == qtr::MolIdType::UID)
                        idToStringWriter << std::make_pair(id, strId);
                    else if (args.molIdType() == qtr::MolIdType::SMILES)
                        idToStringWriter << std::make_pair(id, smiles);
                    else
                        throw std::logic_error("Undefined molecule Id Type");
                    if (propertiesTableWriter != nullptr)
                        *propertiesTableWriter << std::make_pair(id, properties);
                    processed++;
                }
            }
            catch (std::exception &e) {
                LOG(ERROR) << "Fail to parse molecule from " << csvFilePath << " -- " << e.what();
                skipped++;
            }
            if ((processed + skipped) % 100'000 == 0)
                LOG(INFO) << (processed + skipped) << " molecules was processed from " << csvFilePath;
        }
        LOG(INFO) << "Finish processing " << csvFilePath << " : skipped -- " << skipped << ", processed -- "
                  << processed;
        return fingerprintLength;
    }
}

void CSVPreprocessor::run(const PreprocessingArgs &args) {
    ProfileScope("CSV preprocessing");
    std::vector<std::filesystem::path> csvPaths = qtr::findFiles(args.sourceDir(), ".csv");
    std::vector<std::future<uint64_t >> tasks;
    std::atomic_uint64_t counter = 0;
    std::mutex mutex;
    tasks.reserve(csvPaths.size());
    uint64_t fingerprintLength = 0;
    for (auto &csvFile: csvPaths)
        tasks.emplace_back(
                std::async(std::launch::async, parseCSV, std::cref(csvFile), std::cref(args), std::ref(counter),
                           std::ref(mutex)));

    for (auto &task: tasks) {
        uint64_t taskFingerprintLength = task.get();
        if (fingerprintLength == 0)
            fingerprintLength = taskFingerprintLength;
        else
            assert(fingerprintLength == taskFingerprintLength);
    }
    ofstream fingerprintLengthWriter(args.fingerprintLengthFile());
    fingerprintLengthWriter << fingerprintLength;
}
