#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include <glog/logging.h>

#include "Utils.h"
#include "smiles_table_io/SmilesTableReader.h"

#include <map>

ABSL_FLAG(std::string, smiles_dir_path, "",
          "Path to directory with smiles to analyze");

ABSL_FLAG(std::string, stat_file_path, "",
          "Path to file where statistics should be stored");

struct Args {
    std::filesystem::path smilesDirPath;
    std::filesystem::path statFilePath;

    Args(int argc, char *argv[]) {
        absl::ParseCommandLine(argc, argv);
        smilesDirPath = absl::GetFlag(FLAGS_smiles_dir_path);
        qtr::emptyArgument(smilesDirPath, "Please specify smiles_dir_path option");
        statFilePath = absl::GetFlag(FLAGS_stat_file_path);
        qtr::emptyArgument(statFilePath, "Please specify stat_file_path option");
    }
};

void initLogging(int argc, char *argv[]) {
    google::InitGoogleLogging(argv[0]);
    google::SetLogDestination(google::INFO, "analyze_smiles_for_huffman.log");
    FLAGS_alsologtostderr = true;
}

int main(int argc, char *argv[]) {
    initLogging(argc, argv);
    Args args(argc, argv);

    std::map<char, uint64_t> frequency;

    for (auto &smilesTableFile: qtr::findFiles(args.smilesDirPath, qtr::smilesTableExtension)) {
        LOG(INFO) << "Start processing " << smilesTableFile;
        for (const auto &[_, smiles]: qtr::SmilesTableReader(smilesTableFile)) {
            for (char symbol: smiles) {
                frequency[symbol]++;
            }
        }
        LOG(INFO) << "Finish processing " << smilesTableFile;
    }

    std::ofstream out(args.statFilePath);
    for (auto &[symbol, rate]: frequency) {
        out << symbol << ' ' << rate << '\n';
    }

    return 0;
}