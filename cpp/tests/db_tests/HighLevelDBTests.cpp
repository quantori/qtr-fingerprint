#include "gtest/gtest.h"

#include "../utils/TmpDirFixture.h"

#include "PreprocessingArgs.h"
#include "PreprocessorFactory.h"
#include "BuildArgs.h"
#include "DatabaseBuilderFactory.h"
#include "RunArgs.h"
#include "SearchDataLoader.h"
#include "RunModeFactory.h"

using namespace std;
using namespace qtr;

class HighLevelDatabaseTests : public TmpDirFixture {
private:
    static void runFunction(void(*function)(int, char *[]), initializer_list<string> arguments) {
        const static string firstArg = "program";
        vector<char *> args = {const_cast<char *>(firstArg.c_str())};
        for (const string &arg: arguments) {
            args.push_back(const_cast<char *>(arg.c_str()));
        }
        function(int(args.size()), args.data());
    }

public:
    static void preprocess(int argc, char *argv[]) {
        ASSERT_NO_THROW({
                            PreprocessingArgs args(argc, argv);
                            auto preprocessor = PreprocessorFactory::create(args.preprocessingType());
                            preprocessor->run(args);
                        });
    }

    static void buildDB(int argc, char *argv[]) {
        ASSERT_NO_THROW({
                            BuildArgs args(argc, argv);
                            auto builder = DatabaseBuilderFactory::create(args.dbType());
                            builder->build(args);
                        });
    }

    static void runDB(int argc, char *argv[]) {
        ASSERT_NO_THROW({
                            RunArgs args(argc, argv);
                            auto searchData = SearchDataLoader::load(args);
                            unique_ptr<RunMode> mode = RunModeFactory::create(args, searchData);
                            mode->run();
                        });
    }

    static void preprocess(initializer_list<string> args) {
        runFunction(preprocess, args);
    }

    static void buildDB(initializer_list<string> args) {
        runFunction(buildDB, args);
    }

    static void runDB(initializer_list<string> args) {
        runFunction(runDB, args);
    }

    void runQtrDatabase(const filesystem::path &csvPath,
                        const filesystem::path &queriesPath,
                        const vector<size_t> &expectedSummary,
                        size_t drivesCount,
                        size_t parallelizeDepth,
                        size_t treeDepth,
                        size_t runTreads,
                        size_t ansCount,
                        double timeLimit) const {
        filesystem::path csvDir = DataPathManager::getDataDir() / csvPath;
        filesystem::path carbonQueries = DataPathManager::getDataDir() / queriesPath;
        string dbName = "testDB";
        vector<filesystem::path> dbDrives;
        string dbDrivesStr;
        for (size_t i = 0; i < drivesCount; i++) {
            dbDrives.push_back(TmpDirFixture::getTmpDir() / ("drive" + to_string(i)));
            dbDrivesStr += dbDrives.back().string();
            if (i + 1 < drivesCount)
                dbDrivesStr += ',';
        }

        filesystem::path preprocessedDir = TmpDirFixture::getTmpDir() / "preprocessed";
        filesystem::path otherDataDir = TmpDirFixture::getTmpDir() / "other_data";
        filesystem::path summaryFile = TmpDirFixture::getTmpDir() / "summary.txt";

        if (!is_directory(csvDir) || !is_regular_file(carbonQueries)) {
            GTEST_SKIP() << "No data for test";
        }
        for (const auto &dir: {preprocessedDir, otherDataDir, dbDrives[0], dbDrives[1]}) {
            create_directory(dir);
        }
        preprocess({"--preprocessingType=CSV",
                    "--sourceDir=" + csvDir.string(),
                    "--destDir=" + preprocessedDir.string(),
                    "--properties=false",
                    "--molIdType=OriginalId"
                   });


        buildDB({"--dbName=" + dbName,
                 "--dbType=QtrRam",
                 "--sourceDir=" + preprocessedDir.string(),
                 "--destDirs=" + dbDrivesStr,
                 "--otherDestDir=" + otherDataDir.string(),
                 "--parallelizeDepth=" + to_string(parallelizeDepth),
                 "--treeDepth=" + to_string(treeDepth),
                 "--properties=false",
                });

        runDB({"--dbName=" + dbName,
               "--dataDirs=" + dbDrivesStr,
               "--otherDataDir=" + otherDataDir.string(),
               "--threads=" + to_string(runTreads),
               "--mode=FromFile",
               "--ansCount=" + to_string(ansCount),
               "--timeLimit=" + to_string(timeLimit),
               "--queriesFile=" + carbonQueries.string(),
               "--summaryFile=" + summaryFile.string(),
               "--properties=false"
              });

        ifstream summary(summaryFile);
        for (size_t i : expectedSummary) {
            int actual;
            summary >> actual;
            EXPECT_EQ(actual, i);
        }
    }

    void testCarbonDatabase(size_t drivesCount,
                            size_t parallelizeDepth,
                            size_t treeDepth,
                            size_t runTreads,
                            size_t ansCount,
                            double timeLimit) const {
        vector<size_t> summary(200);
        iota(summary.begin(), summary.end(), 1);
        reverse(summary.begin(), summary.end());

        runQtrDatabase("molecules/carbon_csv", "molecules/carbon_queries.txt", summary, drivesCount, parallelizeDepth,
                       treeDepth, runTreads, ansCount, timeLimit);
    }
};

TEST_F(HighLevelDatabaseTests, SimpleCarbonDB) {
    testCarbonDatabase(2, 3, 5, 1, 10000, 2.0);
}

TEST_F(HighLevelDatabaseTests, ParallelCarbonDB) {
    testCarbonDatabase(2, 3, 5, 2, 10000, 2.0);
}
