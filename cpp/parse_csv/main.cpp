#include "absl/flags/flag.h"
#include "absl/flags/parse.h"

#include <filesystem>
#include <future>

#include "Utils.h"
#include "smiles_table_io/SmilesTableWriter.h"
#include "fingerprint_table_io/FingerprintTableWriter.h"
#include "id_to_string_io/IdToStringWriter.h"

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

void parseCSV(const std::filesystem::path &csvFilePath, const Args &args, std::atomic_uint64_t &counter) {
    std::filesystem::path smilesTablePath =
            args.destDirPath / "smilesTables" / (csvFilePath.stem().string() + qtr::smilesTableExtension);
    std::filesystem::path fingerprintTablePath =
            args.destDirPath / "fingerprintTables" / (csvFilePath.stem().string() + qtr::fingerprintTableExtension);
    std::filesystem::path idToStringTablePath =
            args.destDirPath / "idToStringTables" / (csvFilePath.stem().string() + ".csv");

    std::filesystem::create_directory(smilesTablePath.parent_path());
    std::filesystem::create_directory(fingerprintTablePath.parent_path());
    std::filesystem::create_directory(idToStringTablePath.parent_path());

    qtr::SmilesTableWriter smilesTableWriter(smilesTablePath);
    qtr::FingerprintTableWriter fingerprintTableWriter(fingerprintTablePath);
    qtr::IdToStringWriter idToStringWriter(idToStringTablePath);


    std::ifstream in(csvFilePath);
    std::string line;
    uint64_t skipped = 0, processed = 0;
    while (std::getline(in, line)) {
        size_t spaceSepPosReversed = std::find(line.rbegin(), line.rend(), ' ') - line.rbegin();
        size_t tabSepPosReversed = std::find(line.rbegin(), line.rend(), '\t') - line.rbegin();
        size_t sepPosReversed = std::min(spaceSepPosReversed, tabSepPosReversed);
        if (sepPosReversed == line.size()) {
            LOG(WARNING) << "Skip invalid line: " << line;
            continue;
        }
        size_t sepPos = line.size() - 1 - sepPosReversed;
        std::string strId = line.substr(sepPos + 1);

        std::string smiles = line.substr(0, sepPos);

        try {
            qtr::IndigoFingerprint fingerprint = qtr::indigoFingerprintFromSmiles(smiles);

            uint64_t id = counter++;
            smilesTableWriter << std::make_pair(id, smiles);
            fingerprintTableWriter << std::make_pair(id, fingerprint);
            idToStringWriter << std::make_pair(id, strId);
            processed++;
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
    for (auto &csvFile: csvPaths) {
        tasks.emplace_back(
                std::async(std::launch::async, parseCSV, std::cref(csvFile), std::cref(args), std::ref(counter)));
    }
    for (auto &task: tasks) {
        task.get();
    }

    timeTicker.tick("Elapsed time");
    return 0;
}