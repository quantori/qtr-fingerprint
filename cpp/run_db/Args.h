#pragma once

#include "ArgsBase.h"

ABSL_DECLARE_FLAG(std::string, dbName);
ABSL_DECLARE_FLAG(std::string, dbType);
ABSL_DECLARE_FLAG(std::vector<std::string>, dataDirs);
ABSL_DECLARE_FLAG(std::string, otherDataDir);
ABSL_DECLARE_FLAG(uint64_t, threads);
ABSL_DECLARE_FLAG(std::string, mode);
ABSL_DECLARE_FLAG(std::string, queriesFile);
ABSL_DECLARE_FLAG(uint64_t, ansCount);
ABSL_DECLARE_FLAG(double, timeLimit);

namespace qtr {

    class Args : public ArgsBase {
    public:
        enum class Mode {
            BadMode,
            Interactive,
            FromFile,
            Web
        };

    private:
        const static inline std::unordered_map<std::string, Mode> _strToMode = {
                {"interactive", Mode::Interactive},
                {"fromFile",    Mode::FromFile},
                {"web",         Mode::Web}
        };

        const static inline auto strToMode = makeStringToEnumFunction(_strToMode, Mode::BadMode);

    ADD_ARGUMENT_WITH_PARSER(DataBaseType, dbType, DataBaseType::BadType, strToDataBaseType)

    ADD_ARGUMENT(std::string, dbName, "")

    ADD_ARGUMENT_WITH_PARSER(std::vector<std::filesystem::path>, dataDirs, {}, vecStrToVecPath)

    ADD_ARGUMENT(std::filesystem::path, otherDataDir, "")

    ADD_ARGUMENT(uint64_t, threads, -1)

    ADD_ARGUMENT_WITH_PARSER(Mode, mode, Mode::BadMode, strToMode)

    ADD_ARGUMENT(std::filesystem::path, queriesFile, "")

    ADD_ARGUMENT(uint64_t, ansCount, -1)

    ADD_ARGUMENT(double, timeLimit, -1)

    public:
        Args(int argc, char *argv[]) : ArgsBase(argc, argv) {
            parseAndCheck_dbType();
            parseAndCheck_dbName();
            parseAndCheck_dataDirs();
            parseAndCheck_threads();
            parseAndCheck_mode();
            parseAndCheck_ansCount();
            parseAndCheck_timeLimit();

            if (dbType() == ArgsBase::DataBaseType::QtrRam ||
                dbType() == ArgsBase::DataBaseType::QtrDrive) {
                parseAndCheck_otherDataDir();
            }

            if (mode() == Mode::FromFile) {
                parseAndCheck_queriesFile();
            }

            if (dbType() == ArgsBase::DataBaseType::BingoNoSQL) {
                if (mode() != Mode::FromFile) {
                    LOG(ERROR) << "Only FromFile mode is supported for BingoNoSQL";
                    exit(-1);
                }
                if (threads() != 1) {
                    LOG(ERROR) << "Only single-threaded run is supported for BingoNoSQL";
                    exit(-1);
                }
                if (dataDirs().size() != 1) {
                    LOG(ERROR) << "Multiple data folders is not supported for BingoNoSQL";
                    exit(-1);
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
    };
}