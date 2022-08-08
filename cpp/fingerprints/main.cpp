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
#include "RawBucketsIO.h"

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

IndigoFingerprint cutZeroColumns(FullIndigoFingerprint fingerprint) {
    IndigoFingerprint cutFingerprint;
    std::vector<int> zeroColumns = {};
    ifstream fin("zero_columns");
    int number;
    while (fin >> number) {
        zeroColumns.push_back(number);
    }
    int currentZero = 0;
    for (int i = 0; i < fromBytesToBits(fingerprint.sizeInBytes); ++i) {
        if (currentZero < zeroColumns.size() && i == zeroColumns[currentZero]) {
            currentZero++;
            continue;
        }
        cutFingerprint[i - currentZero] = fingerprint[i];
    }
    return cutFingerprint;
}

void createFingerprintCSVFromFile(const string &sdfFile) {
    auto indigoSessionPtr = IndigoSession::create();
    ofstream fout(sdfFile + ".csv");
    for (auto &mol: indigoSessionPtr->iterateSDFile(sdfFile)) {
        try {
            ostringstream fingerprint_line;
            mol->aromatize();
            fingerprint_line << mol->smiles() << " ";
            int fingerprint = indigoFingerprint(mol->id(), "sub");
            hexToBin(indigoToString(fingerprint), fingerprint_line);
            fout << fingerprint_line.str() << endl;
        }
        catch (const exception &e) {
            cerr << e.what();
        }
    }
}

void createRBFromSDF(const filesystem::path &sdfFilePath, const filesystem::path &rbFilePath) {
    auto indigoSessionPtr = IndigoSession::create();
    RawBucketWriter writer(rbFilePath);
    uint64_t skipped = 0;
    for (auto &mol: indigoSessionPtr->iterateSDFile(sdfFilePath)) {
        try {
            mol->aromatize();
        }
        catch (const exception &e) {
            cerr << "Aromatize error: " << e.what() << endl;
            skipped++;
        }
        try {
            int fingerprint = indigoFingerprint(mol->id(), "sub");
            FullIndigoFingerprint fullFingerprints(indigoToString(fingerprint));
            IndigoFingerprint cutFingerprint = cutZeroColumns(fullFingerprints);
            writer.write({mol->smiles(), cutFingerprint});
        }
        catch (const exception &e) {
            cerr << "Other error" << e.what() << endl;
            skipped++;
        }
    }
    cerr << "skipped: " << skipped << endl;
}

ABSL_FLAG(std::string, path_to_sdf_dir, "",
          "Path to dir with sdf files");

ABSL_FLAG(std::string, path_to_rb_dir, "",
          "Path to dir with rb (raw bucket) files");

int main(int argc, char *argv[]) {
    google::InitGoogleLogging(argv[0]);
    absl::ParseCommandLine(argc, argv);
    filesystem::path pathToSdfDir = absl::GetFlag(FLAGS_path_to_sdf_dir);
    filesystem::path pathToRbDir = absl::GetFlag(FLAGS_path_to_rb_dir);
    emptyArgument(pathToSdfDir, "Please specify path_to_sdf_dir option");
    emptyArgument(pathToRbDir, "Please specify path_to_rb_dir option");
    vector<filesystem::path> sdfFiles = findFiles(pathToSdfDir, ".sdf");
    auto startTime = std::chrono::high_resolution_clock::now();
    #pragma omp parallel for
    for (int i = 0; i < sdfFiles.size(); ++i) {
        auto &sdfFilePath = sdfFiles[i];
        auto rbFilePath = pathToRbDir / (string(sdfFilePath.stem()) + ".rb");
        createRBFromSDF(sdfFiles[i], rbFilePath);
    }
    auto endTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = endTime - startTime;
    std::cout << "Elapsed time: " << elapsed_seconds.count() << "s\n";
    return 0;
}
