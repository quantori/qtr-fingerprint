#include "buildBingoNoSQL.h"

#include "bingo-nosql.h"
#include "IndigoMolecule.h"
#include "IndigoSDFileIterator.h"
#include "IndigoException.h"
#include "string_table_io/StringTableReader.h"

#include <future>

using namespace std;
using namespace indigo_cpp;

namespace qtr {
    namespace {
        void storeSmilesTableInDB(int db, const shared_ptr <IndigoSession> &indigoSessionPtr,
                                  const filesystem::path &smilesTablePath) {
            LOG(INFO) << "Start " << smilesTablePath << " parsing";
            IndigoSDFileIterator iterator = indigoSessionPtr->iterateSDFile(smilesTablePath.c_str());
            size_t processedNumber = 0;
            size_t failuresNumber = 0;
            for (const auto &[id, smiles]: StringTableReader(smilesTablePath)) {
                unique_ptr<IndigoMolecule> molecule = nullptr;
                try {
                    molecule = make_unique<IndigoMolecule>(indigoSessionPtr->loadMolecule(smiles));
                    molecule->aromatize();
                    indigoSessionPtr->_checkResult(bingoInsertRecordObj(db, molecule->id()));
                }
                catch (const IndigoException &e) {
                    LOG(ERROR) << "Indigo error: " << e.what();
                    failuresNumber++;
                }
                processedNumber++;
                if (processedNumber % 1000 == 0) {
                    LOG(INFO) << "Processed " << processedNumber << " molecules from " << smilesTablePath;
                }
            }
            LOG(INFO) << "Finish " << smilesTablePath << " parsing " <<
                      "(" << processedNumber << " processed, " << failuresNumber << " failures)";
        }
    }

    void buildBingoNoSQLDb(const Args &args, TimeMeasurer &statisticCollector) {
        TimeMeasurer::FunctionExecutionTimer timer(statisticCollector, "bingo NO SQL building");

        auto indigoSessionPtr = IndigoSession::create();

        auto destDirectories = args.destDirs();
        if (destDirectories.size() > 1) {
            LOG(ERROR) << "Distributed Bingo No SQL Database is not supported. Please, provide only path to " FLAG_NAME(
                    destDirs);
            exit(-1);
        }
        filesystem::path destDir = args.destDirs()[0];
        if (!filesystem::is_directory(destDir)) {
            LOG(ERROR) << "Destination directory (" << destDir << ") does not exist";
            exit(-1);
        }
        int db = bingoCreateDatabaseFile((destDir / args.dbName()).c_str(), "molecule", "");
        vector<future<void>> tasks;
        for (const auto &entry: filesystem::directory_iterator(args.smilesSourceDir())) {
            tasks.push_back(async(launch::async, storeSmilesTableInDB, db, cref(indigoSessionPtr),
                                  entry.path()));
        }
        for (auto &task: tasks) {
            task.wait();
        }
        bingoCloseDatabase(db);
    }
}