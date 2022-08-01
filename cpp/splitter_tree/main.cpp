#include "indigo.h"

#include "IndigoMolecule.h"
#include "IndigoSession.h"
#include "IndigoWriteBuffer.h"
#include "IndigoSDFileIterator.h"

#include <glog/logging.h>

#include <chrono>
#include <iostream>

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"

#include "Utils.h"
#include "Fingerprint.h"
#include "SplitterTree.h"
#include "ColumnsChooser.h"

using namespace indigo_cpp;
using namespace qtr;
using namespace std;

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

/**
 * Merges many sdf files into one, named "0". Suitable for Splitter Tree.
 * @param dir -- directory with sdf files
 */
void prepareSDFsForSplitterTree(const string &dir) {
    auto indigoSessionPtr = IndigoSession::create();
    string filename = dir + "/buckets/0";
    if (dir.back() == '/')
        filename = dir + "buckets/0";
    ofstream fout(filename);
    uint64_t cntMols = 0;
    uint64_t cntSkipped = 0;
    fout.write((char *) (&cntMols), sizeof(cntMols)); // Reserve space for bucket size
    for (auto &sdfFile: findFiles(dir, ".sdf")) {
        for (auto &mol: indigoSessionPtr->iterateSDFile(sdfFile)) {
            try {
                mol->aromatize();
                int fingerprint = indigoFingerprint(mol->id(), "sub");
                FullIndigoFingerprint fp(indigoToString(fingerprint));
                IndigoFingerprint cutFP = cutZeroColumns(fp);
                cutFP.saveBytes(fout);
                fout << mol->smiles() << '\n';
                ++cntMols;
            }
            catch (...) {
                cntSkipped++;
            }
        }
    }
    fout.seekp(0, ios::beg); // Seek to the beginning of file
    fout.write((char *) (&cntMols), sizeof(cntMols));
    if (cntSkipped)
        std::cerr << "Skipped " << cntSkipped << " molecules\n";
}

ABSL_FLAG(std::string, path_to_dir, "",
          "Path to dir with sdf files");

int main(int argc, char *argv[]) {
    google::InitGoogleLogging(argv[0]);
    absl::ParseCommandLine(argc, argv);
    std::string pathToDir = absl::GetFlag(FLAGS_path_to_dir);
    emptyArgument(pathToDir, "Please specify path_to_dir option");
    auto startTime = std::chrono::high_resolution_clock::now();
    prepareSDFsForSplitterTree(pathToDir);
    std::chrono::duration<double> mediumSec = std::chrono::high_resolution_clock::now() - startTime;
    std::cout << "Done preparing in time: " << mediumSec.count() << "s\n";
    SplitterTree tree("/home/Vsevolod.Vaskin/qtr-fingerprint/data/buckets/", "0");
    tree.split(30, 100);
    ofstream fout("/home/Vsevolod.Vaskin/qtr-fingerprint/data/tree");
    tree.saveTo(fout);
    auto chooser = qtr::ColumnsChooser("/home/Vsevolod.Vaskin/qtr-fingerprint/data/buckets/", qtr::correlationColumnsChoose);
    chooser.choose();
    auto endTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = endTime - startTime;
    std::cout << "Elapsed time: " << elapsed_seconds.count() << "s\n";
    return 0;
}
