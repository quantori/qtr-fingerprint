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
};

TEST_F(HighLevelDatabaseTests, CarbonDB) {
    filesystem::path csvDir = DataPathManager::getDataDir() / "molecules/carbon_csv";
    filesystem::path preprocessedDir = TmpDirFixture::getTmpDir() / "preprocessed";
    string dbName = "CarbonDB";
    vector<filesystem::path> dbDrives = {
            TmpDirFixture::getTmpDir() / "drive1",
            TmpDirFixture::getTmpDir() / "drive2"
    };
    filesystem::path otherDataDir = TmpDirFixture::getTmpDir() / "other_data";

    if (!is_directory(csvDir)) {
        GTEST_SKIP() << "No data for test";
    }
    for (const auto& dir : {preprocessedDir, otherDataDir, dbDrives[0], dbDrives[1]}) {
        create_directory(dir);
    }
    preprocess({"--preprocessingType=CSV",
                "--preprocessDir=" + csvDir.string(),
                "--destDir=" + preprocessedDir.string(),
                "--preprocessProperties=false",
               });


    buildDB({"--dbName=" + dbName,
             "--dbType=QtrRam",
             "--sourceDir=" + preprocessedDir.string(),
             "--destDirs=" + dbDrives[0].string() + "," + dbDrives[1].string(),
             "--otherDestDir=" + otherDataDir.string(),
             "--parallelizeDepth=3",
             "--treeDepth=5",
             "--buildProperties=false",
             });
}

