#include "RDKitDatabaseBuilder.h"

#include "Profiling.h"
#include "string_table_io/StringTableReader.h"

#include <future>

using namespace std;

namespace qtr {
    namespace {
        void storeSmilesTableInPickle(std::ofstream &pickleOstream, mutex &mutex,
                                      const filesystem::path &smilesTablePath) {
            LOG(INFO) << "Start " << smilesTablePath << " parsing";
            size_t processedNumber = 0;
            size_t failuresNumber = 0;

            // TODO: maybe use MOLBundle?
            // TODO: use XQMol instead of ROMol?

            for (const auto &[id, smiles]: StringTableReader(smilesTablePath)) {
                try {
                    // TODO uncomment
//                    std::shared_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol(smiles));
                    // TODO: check aromatize
                    string pickle;
                    // TODO: uncomment
                    //  RDKit::MolPickler::pickleMol(*mol, pickle);
                    {
                        lock_guard<std::mutex> guard(mutex);
                        pickleOstream << pickle;
                    }
                    ++processedNumber;
                } catch (const std::exception &e) {
                    LOG(ERROR) << "Indigo error: " << e.what();
                    failuresNumber++;
                }

            }
            pickleOstream.close();
            LOG(INFO) << "Finish " << smilesTablePath << " parsing " <<
                      "(" << processedNumber << " processed, " << failuresNumber << " failures)";
        }
    }

    void RDKitDatabaseBuilder::build(const qtr::BuildArgs &args) {
        ProfileScope("RDKit building");

        auto destDirectories = args.destDirs();
        if (destDirectories.size() > 1) {
            LOG_ERROR_AND_EXIT(
                    "Distributed RDKit Database is not supported. Please, provide only path to " FLAG_NAME(destDirs));
        }
        filesystem::path destDir = args.destDirs()[0];
        if (!filesystem::is_directory(destDir)) {
            LOG_ERROR_AND_EXIT("Destination directory (" + destDir.string() + ") does not exist");
        }
        vector<future<void>> tasks;
        for (const auto &entry: filesystem::directory_iterator(args.smilesSourceDir())) {
            tasks.push_back(async(launch::async, storeSmilesTableInPickle, db, cref(indigoSessionPtr),
                                  entry.path()));
        }
        for (auto &task: tasks) {
            task.wait();
        }
    }
}
