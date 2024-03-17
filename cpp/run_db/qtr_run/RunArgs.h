#pragma once

#include "ArgsBase.h"
#include "modes/RunMode.h"

ABSL_DECLARE_FLAG(std::string, dbName);
ABSL_DECLARE_FLAG(std::string, dbType);
ABSL_DECLARE_FLAG(std::vector<std::string>, dataDirs);
ABSL_DECLARE_FLAG(std::string, otherDataDir);
ABSL_DECLARE_FLAG(uint64_t, threads);
ABSL_DECLARE_FLAG(std::string, mode);
ABSL_DECLARE_FLAG(std::string, queriesFile);
ABSL_DECLARE_FLAG(uint64_t, ansCount);
ABSL_DECLARE_FLAG(double, timeLimit);
ABSL_DECLARE_FLAG(std::string, summaryFile);
ABSL_DECLARE_FLAG(bool, properties);
ABSL_DECLARE_FLAG(bool, verificationStage);
ABSL_DECLARE_FLAG(bool, fingerprintProvided);

namespace qtr {

    class RunArgs : public ArgsBase {
    private:
        const static inline std::unordered_map<std::string, RunMode::Type> _strToMode = {
                {"Interactive", RunMode::Type::Interactive},
                {"FromFile",    RunMode::Type::FromFile},
                {"Web",         RunMode::Type::Web}
        };

        const static inline auto strToMode = makeStringToEnumFunction(_strToMode, RunMode::Type::BadType);

    ADD_ARGUMENT_WITH_PARSER(DatabaseType, dbType, DatabaseType::BadType, strToDataBaseType)

    ADD_ARGUMENT(std::string, dbName, "")

    ADD_ARGUMENT_WITH_PARSER(std::vector<std::filesystem::path>, dataDirs, {}, vecStrToVecPath)

    ADD_ARGUMENT(std::filesystem::path, otherDataDir, "")

    ADD_ARGUMENT(uint64_t, threads, -1)

    ADD_ARGUMENT_WITH_PARSER(RunMode::Type, mode, RunMode::Type::BadType, strToMode)

    ADD_ARGUMENT(std::filesystem::path, queriesFile, "")

    ADD_ARGUMENT(uint64_t, ansCount, -1)

    ADD_ARGUMENT(double, timeLimit, -1)

    ADD_ARGUMENT(std::filesystem::path, summaryFile, "")

    ADD_ARGUMENT(bool, properties, true);

    ADD_ARGUMENT(bool, verificationStage, true);

    ADD_ARGUMENT(bool, fingerprintProvided, false);

    public:
        RunArgs(int argc, char *argv[]) : ArgsBase(argc, argv) {
            parseAndCheck_dbType();
            parseAndCheck_dbName();
            parseAndCheck_dataDirs();
            parseAndCheck_threads();
            parseAndCheck_mode();
            parseAndCheck_ansCount();
            parseAndCheck_timeLimit();
            parse_properties();
            parse_verificationStage();
            parse_fingerprintProvided();

            if (dbType() == DatabaseType::QtrRam ||
                dbType() == DatabaseType::QtrDrive ||
                dbType() == DatabaseType::QtrEnumeration) {
                parseAndCheck_otherDataDir();
            }

            if (dbType() == DatabaseType::QtrEnumeration) {
                if (verificationStage())
                    logErrorAndExit("Verification stage is not implemented for QtrEnumeration");
                if (properties())
                    logErrorAndExit("Properties are not supported for QtrEnumeration");
            }
            if (dbType() == DatabaseType::QtrEnumeration && verificationStage()) {
            }

            if (mode() == RunMode::Type::FromFile) {
                parseAndCheck_queriesFile();
                parse_summaryFile();
            }
            if (mode() != RunMode::Type::FromFile && fingerprintProvided()) {
                logErrorAndExit("fingerprintProvided option is supported for FromFile mode only");
            }

            if (dbType() == DatabaseType::BingoNoSQL) {
                if (properties()) {
                    logErrorAndExit("Only databases without properties are supported for BingoNoSQL");
                }
                if (threads() != 1) {
                    logErrorAndExit("Only single-threaded run is supported for BingoNoSQL");
                }
                if (dataDirs().size() != 1) {
                    logErrorAndExit("Multiple data folders is not supported for BingoNoSQL");
                }
                if (!verificationStage()) {
                    logErrorAndExit("verificationStage=false is not supported for BingoNoSQL");
                }
            }
        }

        [[nodiscard]] inline std::vector<std::filesystem::path> dbDataDirs() const {
            std::vector<std::filesystem::path> res;
            for (auto &dir: dataDirs()) {
                res.emplace_back(dir / dbName());
            }
            return res;
        }

        [[nodiscard]] inline std::filesystem::path dbOtherDataDir() const {
            return otherDataDir() / dbName();
        }

        [[nodiscard]] inline std::filesystem::path ballTreePath() const {
            return dbOtherDataDir() / "tree";
        }

        [[nodiscard]] inline std::filesystem::path smilesTablePath() const {
            return dbOtherDataDir() / "smilesTable";
        }

        [[nodiscard]] inline std::filesystem::path huffmanCoderPath() const {
            return dbOtherDataDir() / "huffman";
        }

        [[nodiscard]] inline std::filesystem::path idToStringDir() const {
            return dbOtherDataDir() / "idToString";
        }

        [[nodiscard]] inline std::filesystem::path propertyTablePath() const {
            return dbOtherDataDir() / "propertyTable";
        }

        [[nodiscard]] inline std::filesystem::path fingerprintLengthFile() const {
            return dbOtherDataDir() / "fingerprintLength";
        }
    };
}
