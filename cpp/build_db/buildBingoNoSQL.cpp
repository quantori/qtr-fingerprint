#include "buildBingoNoSQL.h"

#include "bingo-nosql.h"
#include "IndigoMolecule.h"
#include "IndigoSDFileIterator.h"
#include "IndigoException.h"

#include <future>

using namespace std;
using namespace indigo_cpp;

namespace qtr {
    namespace {
        void storeSDFDataInDB(int db, mutex &lock, const shared_ptr <IndigoSession> &indigoSessionPtr,
                              const filesystem::path &sdfFilePath) {
            LOG(INFO) << "Start " << sdfFilePath << " parsing";
            IndigoSDFileIterator iterator = indigoSessionPtr->iterateSDFile(sdfFilePath.c_str());
            size_t processedNumber = 0;
            size_t failuresNumber = 0;
            for (auto &molecule: iterator) {
                try {
                    molecule->aromatize();
                }
                catch (const IndigoException &e) {
                    LOG(ERROR) << "Aromatize error: " << e.what();
                    failuresNumber++;
                }
                {
                    lock_guard<mutex> lockGuard(lock);
                    bingoInsertRecordObj(db, molecule->id());
                }
                processedNumber++;
                if (processedNumber % 1000 == 0) {
                    LOG(INFO) << "Processed " << processedNumber << " molecules from " << sdfFilePath;
                }
            }

            LOG(INFO) << "Finish " << sdfFilePath << " parsing " <<
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
        mutex indigoLock;
        vector<future<void>> tasks;
        for (const auto &entry: filesystem::directory_iterator(args.sourceDir())) {
            tasks.push_back(async(launch::async, storeSDFDataInDB, db, ref(indigoLock), cref(indigoSessionPtr),
                                  entry.path()));
        }
        for (auto &task: tasks) {
            task.wait();
        }
    }
}