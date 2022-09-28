#include "indigo.h"

#include "IndigoMolecule.h"
#include "IndigoSession.h"
#include "IndigoWriteBuffer.h"
#include "IndigoSDFileIterator.h"

#include <glog/logging.h>

#include <chrono>
#include <unordered_map>
#include <fstream>
#include <iostream>
#include <sstream>

#include <omp.h>
#include "absl/flags/flag.h"
#include "absl/flags/parse.h"

#include "Utils.h"
#include "Fingerprint.h"
#include "raw_bucket_io/RawBucketWriter.h"
#include "columns_io/IndexesReader.h"
#include "columns_io/IndexesWriter.h"

using namespace indigo_cpp;
using namespace qtr;
using namespace std;

static void hexToBin(const char *hexdec, ostringstream &out) {
    const static std::unordered_map<char, std::string> HEX_TO_DEC = {
            {'0', "0000"},
            {'1', "0001"},
            {'2', "0010"},
            {'3', "0011"},
            {'4', "0100"},
            {'5', "0101"},
            {'6', "0110"},
            {'7', "0111"},
            {'7', "0111"},
            {'8', "1000"},
            {'9', "1001"},
            {'a', "1010"},
            {'b', "1011"},
            {'c', "1100"},
            {'d', "1101"},
            {'e', "1110"},
            {'f', "1111"}
    };
    for (int i = 0; hexdec[i]; ++i)
        out << HEX_TO_DEC.at(hexdec[i]);
}

std::vector<size_t> zeroColumns;

IndigoFingerprint cutZeroColumns(FullIndigoFingerprint fingerprint) {
    IndigoFingerprint cutFingerprint;
    int currentZero = 0;
    for (int i = 0; i < fromBytesToBits(FullIndigoFingerprint::sizeInBytes); ++i) {
        if (currentZero < zeroColumns.size() && i == zeroColumns[currentZero]) {
            currentZero++;
            continue;
        }
        cutFingerprint[i - currentZero] = fingerprint[i];
    }
    return cutFingerprint;
}

void createCSVFromSDF(const filesystem::path &sdfFilePath, const filesystem::path &csvFilePath) {
    auto indigoSessionPtr = IndigoSession::create();
    ofstream fout(csvFilePath);
    for (auto &mol: indigoSessionPtr->iterateSDFile(sdfFilePath)) {
        try {
            ostringstream fingerprint_line;
            mol->aromatize();
            fingerprint_line << mol->smiles() << " ";
            int fingerprint = indigoFingerprint(mol->id(), "sub");
            hexToBin(indigoToString(fingerprint), fingerprint_line);
            fout << fingerprint_line.str() << endl;
        }
        catch (const exception &e) {
            LOG(INFO) << e.what();
        }
    }
}

void createRBFromSDF(const filesystem::path &sdfFilePath, const filesystem::path &rbFilePath) {
//    if (filesystem::exists(rbFilePath)) {
//        LOG(INFO) << sdfFilePath << " already exist. It was skipped";
//        return;
//    }
    LOG(INFO) << "start creating " << rbFilePath << " from " << sdfFilePath;
    auto indigoSessionPtr = IndigoSession::create();
    RawBucketWriter writer(rbFilePath);
    uint64_t skipped = 0;
    uint64_t written = 0;
    for (auto &mol: indigoSessionPtr->iterateSDFile(sdfFilePath)) {
        try {
            mol->aromatize();
            int fingerprint = indigoFingerprint(mol->id(), "sub");
            FullIndigoFingerprint fullFingerprints(indigoToString(fingerprint));
            IndigoFingerprint cutFingerprint = cutZeroColumns(fullFingerprints);
            writer << std::make_pair(mol->smiles(), cutFingerprint);
            written++;
        }
        catch (const exception &e) {
            LOG(ERROR) << "Fail to parse molecule from " << sdfFilePath << " -- " << e.what();
            skipped++;
        }
        if ((written + skipped) % 100'000 == 0)
            LOG(INFO) << (written + skipped) << " molecules was processed from " << sdfFilePath;
    }
    LOG(INFO) << "Finish creating " << rbFilePath << " : skipped -- " << skipped << ", written -- " << written;
}

ABSL_FLAG(std::string, path_to_sdf_dir, "",
          "Path to dir with sdf files");

ABSL_FLAG(std::string, path_to_store_dir, "",
          "Path to dir with rb (raw bucket) files");

ABSL_FLAG(std::string, path_to_zero_columns, "",
          "Path to file with empty columns");

ABSL_FLAG(std::string, type_of_output, "",
          "Type of output objects. Possible values: rb, csv");

int main(int argc, char *argv[]) {
    google::InitGoogleLogging(argv[0]);
    google::LogToStderr();
    absl::ParseCommandLine(argc, argv);
    filesystem::path pathToSdfDir = absl::GetFlag(FLAGS_path_to_sdf_dir);
    filesystem::path pathToStoreDir = absl::GetFlag(FLAGS_path_to_store_dir);
    filesystem::path pathToZeroColumns = absl::GetFlag(FLAGS_path_to_zero_columns);
    std::string type_of_output = absl::GetFlag(FLAGS_type_of_output);
    emptyArgument(pathToSdfDir, "Please specify path_to_sdf_dir option");
    emptyArgument(pathToStoreDir, "Please specify path_to_store_dir option");
    emptyArgument(pathToZeroColumns, "Please specify path_to_zero_columns option");
    if (type_of_output != "csv" && type_of_output != "rb") {
        LOG(ERROR) << "Please specify type_of_output option with value \"csv\" or \"rb\"";
        exit(-1);
    }

    IndexesReader columnsReader(pathToZeroColumns);
    columnsReader >> zeroColumns;

    vector<filesystem::path> sdfFiles = findFiles(pathToSdfDir, ".sdf");
    auto startTime = std::chrono::high_resolution_clock::now();
#pragma omp parallel for
    for (int i = 0; i < sdfFiles.size(); ++i) {
        auto &sdfFilePath = sdfFiles[i];
        auto outFilePath = pathToStoreDir / (string(sdfFilePath.stem()) + "." + type_of_output);
        if (type_of_output == "csv") {
            createCSVFromSDF(sdfFiles[i], outFilePath);
        } else if (type_of_output == "rb") {
            createRBFromSDF(sdfFiles[i], outFilePath);
        }
    }
    auto endTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = endTime - startTime;
    std::cout << "Elapsed time: " << elapsed_seconds.count() << "s\n";
    return 0;
}
