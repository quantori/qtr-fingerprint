#include <filesystem>
#include <string>
#include <numeric>

#include <glog/logging.h>
#include <absl/flags/flag.h>
#include <absl/flags/parse.h>

#include "Utils.h"
#include "io/BufferedReader.h"
#include "BallTreeRAMSearchEngine.h"
#include "BallTreeDriveSearchEngine.h"
#include "fingerprint_table_io/FingerprintTableReader.h"
#include "IndigoSubstructureMatcher.h"
#include "IndigoQueryMolecule.h"
#include "smiles_table_io/SmilesTableReader.h"

ABSL_FLAG(std::vector<std::string>, data_dir_paths, {},
          "Path to directories where data are stored");

ABSL_FLAG(std::string, other_data_path, "",
          "Path to directory with other files");

ABSL_FLAG(std::string, db_name, "",
          "Name of folders with data base's files");

ABSL_FLAG(uint64_t, start_search_depth, -1,
          "Depth for start search from. There will be 2^start_search_depth threads.");

ABSL_FLAG(std::string, mode, "",
          "possible modes: interactive, from_file");

ABSL_FLAG(std::string, input_file, "",
          "file to load test molecules from");

ABSL_FLAG(std::uint64_t, ans_count, -1,
          "how many answers program will find");

void initLogging(int argc, char *argv[]) {
    google::InitGoogleLogging(argv[0]);
    google::SetLogDestination(google::INFO, "run_db.info");
    FLAGS_alsologtostderr = true;
}

struct Args {

    enum class Mode {
        Interactive,
        FromFile
    };

    std::vector<std::filesystem::path> dataDirPaths;
    std::filesystem::path otherDataPath;
    std::string dbName;
    uint64_t startSearchDepth;
    Mode mode;
    std::filesystem::path inputFile;
    uint64_t ansCount;

    std::vector<std::filesystem::path> dbDataDirsPaths;
    std::filesystem::path dbOtherDataPath;

    std::filesystem::path ballTreePath;
    std::filesystem::path smilesTablePath;
    std::filesystem::path smilesRandomAccessTablePath;
    std::filesystem::path fingerprintTablesPath;

    Args(int argc, char *argv[]) {
        absl::ParseCommandLine(argc, argv);

        std::vector<std::string> dataDirPathsStrings = absl::GetFlag(FLAGS_data_dir_paths);
        std::copy(dataDirPathsStrings.begin(), dataDirPathsStrings.end(), std::back_inserter(dataDirPaths));
        qtr::emptyArgument(dataDirPaths, "Please specify data_dir_paths option");
        for (size_t i = 0; i < dataDirPaths.size(); i++) {
            LOG(INFO) << "dataDirPaths[" << i << "]: " << dataDirPaths[i];
        }

        otherDataPath = absl::GetFlag(FLAGS_other_data_path);
        qtr::emptyArgument(otherDataPath, "Please specify other_data_path option");
        LOG(INFO) << "otherDataPath: " << otherDataPath;

        dbName = absl::GetFlag(FLAGS_db_name);
        qtr::emptyArgument(dbName, "Please specify db_name option");
        LOG(INFO) << "dbName: " << dbName;

        startSearchDepth = absl::GetFlag(FLAGS_start_search_depth);
        if (startSearchDepth == -1) {
            LOG(INFO) << "Please specify start_search_depth option";
            exit(-1);
        }
        LOG(INFO) << "startSearchDepth: " << startSearchDepth;

        std::string modeStr = absl::GetFlag(FLAGS_mode);
        qtr::emptyArgument(modeStr, "Please specify mode option");
        inputFile = absl::GetFlag(FLAGS_input_file);
        LOG(INFO) << "inputFile: " << inputFile;
        if (modeStr == "interactive") {
            mode = Mode::Interactive;
            LOG(INFO) << "mode: interactive";
        } else if (modeStr == "from_file") {
            mode = Mode::FromFile;
            LOG(INFO) << "mode: fromFile";
            qtr::emptyArgument(inputFile, "Please specify input_file option");
        } else {
            LOG(ERROR) << "Bad mode option value";
            exit(-1);
        }

        ansCount = absl::GetFlag(FLAGS_ans_count);
        LOG(INFO) << "ansCount: " << ansCount;

        for (auto &dir: dataDirPaths) {
            dbDataDirsPaths.emplace_back(dir / dbName);
        }
        for (size_t i = 0; i < dbDataDirsPaths.size(); i++) {
            LOG(INFO) << "dbDataDirPaths[" << i << "]: " << dbDataDirsPaths[i];
        }

        dbOtherDataPath = otherDataPath / dbName;
        LOG(INFO) << "dbOtherDataPath" << dbOtherDataPath;

        ballTreePath = dbOtherDataPath / "tree";
        LOG(INFO) << "ballTreePath: " << ballTreePath;

        smilesTablePath = dbOtherDataPath / "smilesTable";
        LOG(INFO) << "smilesTablePath: " << smilesTablePath;

        smilesRandomAccessTablePath = dbOtherDataPath / "smilesRandomAccessTablePath";
        LOG(INFO) << "smilesRandomAccessTablePath: " << smilesRandomAccessTablePath;

        fingerprintTablesPath = dbOtherDataPath / "fingerprintTables";
        LOG(INFO) << "fingerprintTablesPath " << fingerprintTablesPath;
    }
};

template<typename SmilesTable>
std::vector<std::string>
getSmiles(const std::vector<size_t> &smilesIndexes, SmilesTable &smilesTable, uint64_t ansCount) {
    std::vector<std::string> result;
    for (size_t i = 0; i < smilesIndexes.size() && i < ansCount; i++) {
        result.emplace_back(smilesTable[smilesIndexes[i]]);
    }
    return result;
}

void printSmiles(const std::vector<std::string> &smiles) {
    size_t answersToPrint = std::min(size_t(10), smiles.size());
    for (size_t i = 0; i < answersToPrint; i++) {
        if (i == 0)
            std::cout << "[";
        std::cout << '\"' << smiles[i] << '\"';
        std::cout << (i + 1 == answersToPrint ? "]\n" : ", ");
    }
}

template<typename SmilesTable>
bool doSearch(const std::string &querySmiles, const qtr::BallTreeSearchEngine &ballTree,
              SmilesTable &smilesTable, const Args &args) {
    qtr::IndigoFingerprint fingerprint;
    try {
        fingerprint = qtr::IndigoFingerprintFromSmiles(querySmiles);
    }
    catch (std::exception &exception) {
        std::cout << "skip query:" << exception.what() << std::endl;
        return false;
    }

    auto filter = [&smilesTable, &querySmiles](size_t ansId) {
        auto ansSmiles = smilesTable[ansId];
        auto indigoSessionPtr = indigo_cpp::IndigoSession::create();
        auto queryMol = indigoSessionPtr->loadQueryMolecule(querySmiles);
        queryMol.aromatize();
        try {
            auto candidateMol = indigoSessionPtr->loadMolecule(ansSmiles);
            candidateMol.aromatize();
            auto matcher = indigoSessionPtr->substructureMatcher(candidateMol);
            return bool(indigoMatch(matcher.id(), queryMol.id()));
        }
        catch (std::exception &e) {
            LOG(ERROR) << "Error while filtering answer. "
                          "Query: " << querySmiles << ", candidate: " << ansId << " " << ansSmiles << ", error: "
                       << e.what();
            return false;
        }
    };
    auto candidateIndexes = ballTree.search(fingerprint, args.ansCount, args.startSearchDepth, filter);
    auto answerSmiles = getSmiles(candidateIndexes, smilesTable, args.ansCount);
    LOG(INFO) << "found answers: " << answerSmiles.size();
    printSmiles(answerSmiles);
    return true;
}


template<typename SmilesTable>
void runInteractive(const qtr::BallTreeSearchEngine &ballTree, SmilesTable &smilesTable,
                    qtr::TimeTicker &timeTicker, const Args &args) {
    while (true) {
        std::cout << "Enter smiles: ";
        std::string smiles;
        std::cin >> smiles;
        if (smiles.empty())
            break;
        timeTicker.tick();
        if (!doSearch(smiles, ballTree, smilesTable, args))
            continue;
        std::cout << "Search time: " << timeTicker.tick("Search time") << std::endl;
    }
}

template<typename SmilesTable>
void runFromFile(const qtr::BallTreeSearchEngine &ballTree, SmilesTable &smilesTable,
                 qtr::TimeTicker &timeTicker, const Args &args) {
    std::ifstream input(args.inputFile);
    std::vector<std::string> queries;
    while (input.peek() != EOF) {
        std::string query;
        input >> query;
        std::string otherInfoInLine;
        std::getline(input, otherInfoInLine);
        queries.emplace_back(query);
    }
    std::vector<double> times;
    size_t skipped = 0;
    LOG(INFO) << "Loaded " << queries.size() << " queries";
    for (size_t i = 0; i < queries.size(); i++) {
        LOG(INFO) << "Start search for " << i << ": " << queries[i];
        timeTicker.tick();
        if (!doSearch(queries[i], ballTree, smilesTable, args)) {
            skipped++;
            continue;
        }
        times.emplace_back(timeTicker.tick("search molecule " + std::to_string(i) + ": " + queries[i]));
    }

    double mean = std::accumulate(times.begin(), times.end(), 0.0) / times.size();
    double min = *std::min_element(times.begin(), times.end());
    double max = *std::max_element(times.begin(), times.end());
    std::sort(times.begin(), times.end());
    double median = times[times.size() / 2];
    double p60 = times[times.size() * 6 / 10];
    double p70 = times[times.size() * 7 / 10];
    double p80 = times[times.size() * 8 / 10];
    double p90 = times[times.size() * 9 / 10];
    double p95 = times[times.size() * 95 / 100];
    LOG(INFO) << "skipped queries: " << skipped;
    LOG(INFO) << "  mean: " << mean;
    LOG(INFO) << "   max: " << max;
    LOG(INFO) << "   min: " << min;
    LOG(INFO) << "median: " << median;
    LOG(INFO) << "60%: " << p60 << " | 70%: " << p70 << " | 80%: " << p80 << " | 90%: " << p90 << " | 95%: " << p95;
}

void loadSmilesTable(std::vector<std::string> &smilesTable, const std::filesystem::path &smilesTablePath) {
    LOG(INFO) << "Start smiles table loading";
    for (const auto &[id, smiles]: qtr::SmilesTableReader(smilesTablePath)) {
        assert(id == smilesTable.size());
        smilesTable.emplace_back(smiles);
    }
    LOG(INFO) << "Finish smiles table loading";
}


int main(int argc, char *argv[]) {
    initLogging(argc, argv);
    Args args(argc, argv);

    qtr::TimeTicker timeTicker;
    qtr::BufferedReader ballTreeReader(args.ballTreePath);

    std::vector<std::string> smilesTable;
    auto loadSmilesTableTask = std::async(std::launch::async, loadSmilesTable, std::ref(smilesTable),
                                          std::cref(args.smilesTablePath));
    LOG(INFO) << "Start ball tree loading";
    qtr::BallTreeRAMSearchEngine ballTree(ballTreeReader, args.dbDataDirsPaths);
    LOG(INFO) << "Finish ball tree loading";
    loadSmilesTableTask.get();
    timeTicker.tick("DB initialization");
    if (args.mode == Args::Mode::Interactive)
        runInteractive(ballTree, smilesTable, timeTicker, args);
    else if (args.mode == Args::Mode::FromFile)
        runFromFile(ballTree, smilesTable, timeTicker, args);

    return 0;
}