#include "indigo.h"

#include "IndigoMolecule.h"
#include "IndigoSession.h"
#include "IndigoWriteBuffer.h"
#include "IndigoSDFileIterator.h"

#include <glog/logging.h>

#include <chrono>
#include <unordered_map>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

#include <omp.h>
#include "absl/flags/flag.h"
#include "absl/flags/parse.h"

#include "Utils.h"
#include "Fingerprint.h"

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

IndigoFingerprint cutZeroColumns(Fingerprint<467> fingerprint) {
    IndigoFingerprint cutFingerprint;
    std::vector<int> zeroColumns = {};
    ifstream fin("zero_columns");
    int number;
    while (fin >> number) {
        zeroColumns.push_back(number);
    }
    int currentZero = 0;
    for (int i = 0; i < 467 * CHAR_BIT; ++i) {
        if (currentZero < zeroColumns.size() && i == zeroColumns[currentZero]) {
            currentZero++;
            continue;
        }
        cutFingerprint[i - currentZero] = fingerprint[i];
    }
    return cutFingerprint;
}

vector<string> findFiles(const string &pathToDir, const string &extension) {
    vector<string> sdfFiles;
    for (const auto &entry: std::filesystem::recursive_directory_iterator(pathToDir)) {
        if (entry.path().extension() == extension) {
            cout << entry.path().string() << endl;
            sdfFiles.push_back(entry.path().string());
        }
    }
    return sdfFiles;
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

ABSL_FLAG(std::string, path_to_dir, "",
          "Path to dir with sdf files");

int main(int argc, char *argv[]) {
    google::InitGoogleLogging(argv[0]);
    absl::ParseCommandLine(argc, argv);
    std::string pathToDir = absl::GetFlag(FLAGS_path_to_dir);
    emptyArgument(pathToDir, "Please specify path_to_dir option");
    vector<string> sdfFiles = findFiles(pathToDir, ".sdf");
    auto startTime = std::chrono::high_resolution_clock::now();
#pragma omp parallel for
    for (int i = 0; i < sdfFiles.size(); ++i) {
        createFingerprintCSVFromFile(sdfFiles[i]);
    }
    auto endTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = endTime - startTime;
    std::cout << "Elapsed time: " << elapsed_seconds.count() << "s\n";
    return 0;
}
