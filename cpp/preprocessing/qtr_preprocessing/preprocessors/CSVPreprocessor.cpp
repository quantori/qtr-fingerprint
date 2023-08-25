#include <future>
#include "CSVPreprocessor.h"

#include "Utils.h"
#include "string_table_io/StringTableWriter.h"
#include "fingerprint_table_io/FingerprintTableWriter.h"
#include "id_to_string_io/IdToStringWriter.h"
#include "properties_table_io/PropertiesTableWriter.h"
#include "PropertiesFilter.h"
#include "Profiling.h"

namespace qtr {
    namespace {
        std::vector<std::string> splitString(const std::string &str) {
            static const char delimiter = '\t';
            std::vector<std::string> substrings;
            std::size_t start = 0;
            std::size_t end = str.find(delimiter);
            while (end != std::string::npos) {
                substrings.push_back(str.substr(start, end - start));
                start = end + 1;
                end = str.find(delimiter, start);
            }
            substrings.push_back(str.substr(start, end));
            return substrings;
        }

        void
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
            std::filesystem::create_directory(propertyTablesPath.parent_path());

            StringTableWriter smilesTableWriter(smilesTablePath);
            FingerprintTableWriter fingerprintTableWriter(fingerprintTablePath);
            IdToStringWriter idToStringWriter(idToStringTablePath);
            PropertiesTableWriter propertiesTableWriter(propertyTablesPath);

            std::ifstream in(csvFilePath);
            std::string line;
            uint64_t skipped = 0, processed = 0;
            while (std::getline(in, line)) {
                auto lineElements = splitString(line);
                if (lineElements.size() != 2 + PropertiesFilter::Properties().size()) {
                    LOG(WARNING) << "Skip line with invalid (" << lineElements.size()
                                 << ") number of arguments elements: "
                                 << line;
                    skipped++;
                    continue;
                }

                std::string strId = lineElements[0];
                std::string smiles = lineElements[1];
                PropertiesFilter::Properties properties{};
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

                try {
                    IndigoFingerprint fingerprint = indigoFingerprintFromSmiles(smiles);
                    {
                        std::lock_guard lock(mutex);
                        uint64_t id = counter++;
                        smilesTableWriter << std::make_pair(id, smiles);
                        fingerprintTableWriter << std::make_pair(id, fingerprint);
                        idToStringWriter << std::make_pair(id, strId);
                        propertiesTableWriter << std::make_pair(id, properties);
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
        }
    }

    void CSVPreprocessor::run(const PreprocessingArgs &args) {
        ProfileScope("CSV preprocessing");
        std::vector<std::filesystem::path> csvPaths = qtr::findFiles(args.sourceDir(), ".csv");
        std::vector<std::future<void>> tasks;
        std::atomic_uint64_t counter = 0;
        std::mutex mutex;
        for (auto &csvFile: csvPaths) {
            tasks.emplace_back(
                    std::async(std::launch::async, parseCSV, std::cref(csvFile), std::cref(args), std::ref(counter),
                               std::ref(mutex)));
        }
        for (auto &task: tasks) {
            task.get();
        }
    }
} // qtr