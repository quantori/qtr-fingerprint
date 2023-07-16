#pragma once

#include "ArgsBase.h"
#include <absl/flags/declare.h>

ABSL_DECLARE_FLAG(std::string, dbType);

ABSL_DECLARE_FLAG(std::string, dbName);

ABSL_DECLARE_FLAG(std::string, sourceDir);

ABSL_DECLARE_FLAG(std::vector<std::string>, destDirs);

ABSL_DECLARE_FLAG(std::string, otherDestDir);

ABSL_DECLARE_FLAG(uint64_t, parallelizeDepth);

ABSL_DECLARE_FLAG(uint64_t, treeDepth);


namespace qtr {
    class Args : public ArgsBase {
    ADD_ARGUMENT_WITH_PARSER(DatabaseType, dbType, DatabaseType::BadType, strToDataBaseType)

    ADD_ARGUMENT(std::string, dbName, std::string())

    ADD_ARGUMENT(std::filesystem::path, sourceDir, "")

    ADD_ARGUMENT_WITH_PARSER(std::vector<std::filesystem::path>, destDirs, {}, vecStrToVecPath)

    ADD_ARGUMENT(std::filesystem::path, otherDestDir, "")

    ADD_ARGUMENT(uint64_t, parallelizeDepth, 0)

    ADD_ARGUMENT(uint64_t, treeDepth, 0)

    public:
        Args(int argc, char *argv[]) : ArgsBase(argc, argv) {
            parseAndCheck_dbName();
            parseAndCheck_dbType();
            parseAndCheck_sourceDir();
            parseAndCheck_destDirs();
            if (dbType() == DatabaseType::QtrDrive || dbType() == DatabaseType::QtrRam) {
                parseAndCheck_otherDestDir();
                parseAndCheck_parallelizeDepth();
                parseAndCheck_treeDepth();
            } else if (dbType() == DatabaseType::BingoNoSQL) {
                // No specific arguments for BingoNoSQL
            } else {
                LOG(ERROR) << "A case that should not have been executed has been executed";
                exit(-1);
            }
        }

        [[nodiscard]] inline std::filesystem::path smilesSourceDir() const {
            return sourceDir() / "smilesTables";
        }

        [[nodiscard]] inline std::filesystem::path fingerprintTablesSourceDir() const {
            return sourceDir() / "fingerprintTables";
        }

        [[nodiscard]] inline std::filesystem::path idToStringSourceDir() const {
            return sourceDir() / "idToStringTables";
        }

        [[nodiscard]] inline std::filesystem::path propertyTablesSourceDir() const {
            return sourceDir() / "propertyTables";
        }

        [[nodiscard]] inline std::filesystem::path dbOtherDataPath() const {
            return otherDestDir() / dbName();
        }

        [[nodiscard]] inline std::filesystem::path ballTreePath() const {
            return dbOtherDataPath() / "tree";
        }

        [[nodiscard]] inline std::filesystem::path smilesTablePath() const {
            return dbOtherDataPath() / "smilesTable";
        }

        [[nodiscard]] inline std::filesystem::path huffmanCoderPath() const {
            return dbOtherDataPath() / "huffman";
        }

        [[nodiscard]] inline std::filesystem::path idToStringDestinationDir() const {
            return dbOtherDataPath() / "idToString";
        }

        [[nodiscard]] inline std::filesystem::path propertyTableDestinationPath() const {
            return dbOtherDataPath() / "propertyTable";
        }

        [[nodiscard]] inline std::vector<std::filesystem::path> dbDataDirs() const {
            std::vector<std::filesystem::path> res;
            for (auto &dir: destDirs()) {
                res.emplace_back(dir / dbName());
            }
            return res;
        }
    };
}
