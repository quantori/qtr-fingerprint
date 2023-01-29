#include "absl/flags/flag.h"
#include "absl/flags/parse.h"

#include <filesystem>
#include <future>

#include "Utils.h"
#include "smiles_table_io/SmilesTableWriter.h"
#include "fingerprint_table_io/FingerprintTableWriter.h"
#include "id_to_string_io/IdToStringWriter.h"
#include "properties_table_io/PropertiesTableWriter.h"
#include "PropertiesFilter.h"

ABSL_FLAG(std::string, source_dir_path, "",
          "Path to dir with csv files");

ABSL_FLAG(std::string, dest_dir_path, "",
          "Path to directory where parsed data should be stored");

struct Args {
    std::filesystem::path sourceDirPath;
    std::filesystem::path destDirPath;

    Args(int argc, char *argv[]) {
        absl::ParseCommandLine(argc, argv);

        sourceDirPath = absl::GetFlag(FLAGS_source_dir_path);
        qtr::emptyArgument(sourceDirPath, "Please specify source_dir_path option");

        destDirPath = absl::GetFlag(FLAGS_dest_dir_path);
        qtr::emptyArgument(destDirPath, "Please specify dest_dir_path option");
    }
};

std::vector<std::string> splitString(const std::string &str, char delimiter) {
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
parseCSV(const std::filesystem::path &csvFilePath, const Args &args, std::atomic_uint64_t &counter, std::mutex &mutex) {
    std::filesystem::path smilesTablePath =
            args.destDirPath / "smilesTables" / (csvFilePath.stem().string() + qtr::smilesTableExtension);
    std::filesystem::path fingerprintTablePath =
            args.destDirPath / "fingerprintTables" / (csvFilePath.stem().string() + qtr::fingerprintTableExtension);
    std::filesystem::path idToStringTablePath =
            args.destDirPath / "idToStringTables" / (csvFilePath.stem().string() + ".csv");
    std::filesystem::path propertiesTablePath =
            args.destDirPath / "propertiesTable" / csvFilePath.stem();

    std::filesystem::create_directory(smilesTablePath.parent_path());
    std::filesystem::create_directory(fingerprintTablePath.parent_path());
    std::filesystem::create_directory(idToStringTablePath.parent_path());
    std::filesystem::create_directory(propertiesTablePath.parent_path());

    qtr::SmilesTableWriter smilesTableWriter(smilesTablePath);
    qtr::FingerprintTableWriter fingerprintTableWriter(fingerprintTablePath);
    qtr::IdToStringWriter idToStringWriter(idToStringTablePath);
    qtr::PropertiesTableWriter propertiesTableWriter(propertiesTablePath);

    std::ifstream in(csvFilePath);
    std::string line;
    uint64_t skipped = 0, processed = 0;
    while (std::getline(in, line)) {
        auto lineElements = splitString(line, '\t');
        if (lineElements.size() != 2 + qtr::PropertiesFilter::Properties::size()) {
            LOG(WARNING) << "Skip line with invalid (" << lineElements.size() << ") number of arguments elements: "
                         << line;
            skipped++;
            continue;
        }

        std::string strId = lineElements[0];
        std::string smiles = lineElements[1];
        qtr::PropertiesFilter::Properties properties{};
        try {
            for (size_t i = 2, j = 0; j < qtr::PropertiesFilter::Properties::size(); i++, j++) {
                properties[i] = std::stof(lineElements[i]);
            }
        }
        catch (const std::invalid_argument &e) {
            LOG(WARNING) << "Skip line with invalid property: " << line;
            skipped++;
            continue;
        }

        try {
            qtr::IndigoFingerprint fingerprint = qtr::indigoFingerprintFromSmiles(smiles);
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
    LOG(INFO) << "Finish processing " << csvFilePath << " : skipped -- " << skipped << ", processed -- " << processed;
}

int main(int argc, char *argv[]) {
    google::InitGoogleLogging(argv[0]);
    google::LogToStderr();
    Args args(argc, argv);

    qtr::TimeTicker timeTicker;
    std::vector<std::filesystem::path> csvPaths = qtr::findFiles(args.sourceDirPath, ".csv");
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

    timeTicker.tick("Elapsed time");
    return 0;
}