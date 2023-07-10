#pragma once

#include "ArgsBase.h"
#include <absl/flags/declare.h>

ABSL_DECLARE_FLAG(std::string, dbType);

ABSL_DECLARE_FLAG(std::string, dbName);

ABSL_DECLARE_FLAG(std::string, sourceDirPath);

ABSL_DECLARE_FLAG(std::vector<std::string>, destDirPaths);

ABSL_DECLARE_FLAG(std::string, otherDestDirPath);

ABSL_DECLARE_FLAG(uint64_t, parallelizeDepth);

ABSL_DECLARE_FLAG(uint64_t, treeDepth);

namespace qtr {
    class Args : public ArgsBase {
    ADD_ARGUMENT_WITH_PARSER(DataBaseType, dbType, DataBaseType::BadType, strToDataBaseType)

    ADD_ARGUMENT(std::string, dbName, std::string())

    ADD_ARGUMENT(std::filesystem::path, sourceDirPath, "")

    ADD_ARGUMENT_WITH_PARSER(std::vector<std::filesystem::path>, destDirPaths, {}, vecStrToVecPath)

    ADD_ARGUMENT(std::filesystem::path, otherDestDirPath, "")

    ADD_ARGUMENT(uint64_t, parallelizeDepth, 0)

    ADD_ARGUMENT(uint64_t, treeDepth, 0)

    public:
        Args(int argc, char *argv[]) : ArgsBase(argc, argv) {
            parseAndCheck_dbName();
            parseAndCheck_dbType();
            parseAndCheck_sourceDirPath();
            parseAndCheck_destDirPaths();
            if (dbType() == DataBaseType::QtrDrive || dbType() == DataBaseType::QtrRam) {
                parseAndCheck_otherDestDirPath();
                parseAndCheck_parallelizeDepth();
                parseAndCheck_treeDepth();
            } else if (dbType() == qtr::ArgsBase::DataBaseType::BingoNoSQL) {
                // TODO: implement Bingo arguments
            } else {
                LOG(ERROR) << "A case that should not have been executed has been executed";
                exit(-1);
            }
        }

        [[nodiscard]] inline std::filesystem::path smilesSourceDirPath() const {
            return sourceDirPath() / "smilesTables";
        }

        [[nodiscard]] inline std::filesystem::path fingerprintTablesSourceDirPath() const {
            return sourceDirPath() / "fingerprintTables";
        }

        [[nodiscard]] inline std::filesystem::path idToStringSourceDirPath() const {
            return sourceDirPath() / "idTostd::stringTables";
        }

        [[nodiscard]] inline std::filesystem::path propertyTablesSourceDirPath() const {
            return sourceDirPath() / "propertyTables";
        }

        [[nodiscard]] inline std::filesystem::path dbOtherDataPath() const {
            return otherDestDirPath() / dbName();
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

        [[nodiscard]] inline std::filesystem::path idToStringDestinationDirPath() const {
            return dbOtherDataPath() / "idTostd::string";
        }

        [[nodiscard]] inline std::filesystem::path propertyTableDestinationPath() const {
            return dbOtherDataPath() / "propertyTable";
        }

        [[nodiscard]] inline std::vector<std::filesystem::path> dbDataDirPaths() const {
            std::vector<std::filesystem::path> res;
            for (auto &dir: destDirPaths()) {
                res.emplace_back(dir / dbName());
            }
            return res;
        }
    };
}
