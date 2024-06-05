#include "RDKitDatabaseBuilder.h"

#include "Profiling.h"
#include "string_table_io/StringTableReader.h"

#include "GraphMol/GraphMol.h"
#include "GraphMol/SmilesParse/SmilesParse.h"
#include "GraphMol/MolPickler.h"
#include "GraphMol/SubstructLibrary/SubstructLibrary.h"

#include <future>

using namespace std;

namespace qtr {
    namespace {
        void storeSmilesTableInPickle(const filesystem::path &moleculesDir,
                                      const filesystem::path &smilesTablePath) {
            LOG(INFO) << "Start " << smilesTablePath << " parsing";
            size_t processedNumber = 0;
            size_t failuresNumber = 0;

            std::ofstream out(moleculesDir / (smilesTablePath.stem().string() + ".pkl"), ios::binary);
            for (const auto &[id, smiles]: StringTableReader(smilesTablePath)) {
                try {
                    std::shared_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol(smiles));
                    if (mol == nullptr) {
                        throw std::runtime_error("Cannot parse molecule: " + smiles);
                    }
                    RDKit::MolPickler::pickleMol(*mol, out);
                    ++processedNumber;
                } catch (const std::exception &e) {
                    LOG(ERROR) << "Indigo error: " << e.what();
                    failuresNumber++;
                }

            }
            LOG(INFO) << "Finish " << smilesTablePath << " parsing " <<
                      "(" << processedNumber << " processed, " << failuresNumber << " failures)";
        }
    }

    void RDKitDatabaseBuilder::build(const qtr::BuildArgs &args) {
        ProfileScope("RDKit building");

        auto destDirectories = args.dbDataDirs();
        if (destDirectories.size() > 1) {
            LOG_ERROR_AND_EXIT(
                    "Distributed RDKit Database is not supported. Please, provide only path to " FLAG_NAME(destDirs));
        }
        filesystem::path destDir = args.dbDataDirs()[0];
        filesystem::create_directory(destDir);
        if (!filesystem::is_directory(destDir)) {
            LOG_ERROR_AND_EXIT("Destination directory (" + destDir.string() + ") does not exist");
        }
        vector<future<void>> tasks;
        filesystem::path moleculesDir = destDir / "molecules";
        filesystem::create_directory(moleculesDir);
        for (const auto &entry: filesystem::directory_iterator(args.smilesSourceDir())) {
            tasks.push_back(
                    async(launch::async, storeSmilesTableInPickle, cref(moleculesDir), entry.path()));
        }

        auto fingerprintsDir = destDir / "fingerprintTables";
        filesystem::create_directory(fingerprintsDir);
        for (const auto &entry: filesystem::directory_iterator(args.fingerprintTablesSourceDir())) {
            tasks.push_back(async(launch::async, [entry, fingerprintsDir] {
                filesystem::copy_file(entry.path(), fingerprintsDir / entry.path().filename());
            }));
        }
        for (auto &task: tasks) {
            task.wait();
        }
    }
}
