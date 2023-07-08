#pragma once

#include <glog/logging.h>
#include <absl/flags/flag.h>
#include <absl/flags/parse.h>

#include <string>
#include <vector>
#include <filesystem>
#include <functional>
#include <map>

// TODO: сделать универсальную систему для работы с аргументами

#define OPTION_NAME(x) #x
#define GET_FLAG(x) absl::GetFlag(FLAGS_##x)

#define ADD_ARGUMENT_WITH_PARSER(type, name, emptyVal, parseFunction) \
    protected: \
        type _##name = emptyVal; \
        inline void check_##name() const { \
            if (_##name == emptyVal) { \
                LOG(INFO) << "Please specify " << #name << " option"; \
                exit(0); \
            } \
        } \
        inline void parse_##name() {       \
            _##name = parseFunction(absl::GetFlag(FLAGS_##name)); \
            \
        } \
    public: \
        [[nodiscard]] inline type name() const { \
            if (_##name == emptyVal) { \
                LOG(ERROR) << "Try to get uninitialized  field " << #name; \
                exit(-1); \
            } \
            return _##name; \
        }



#define ADD_ARGUMENT(type, name, emptyVal) ADD_ARGUMENT_WITH_PARSER(type, name, emptyVal, std::identity())

ABSL_FLAG(int, test, 0,
          "Name of folders with data base's files");

namespace qtr {

    class ArgsBase {
    public:
        enum class DataBaseType {
            BadType,
            QtrDrive,
            QtrRam,
            BingoNoSQL // TODO implement
        };

        static inline const std::unordered_map<std::string, DataBaseType> strToDataBaseType =
                {{"QtrDrive", DataBaseType::QtrDrive},
                 {"QrtRam", DataBaseType::QtrRam},
                 {"BingoNoSQL", DataBaseType::BingoNoSQL}};

        // todo: перенести map в protected зону, завершить функцию парсинга databasetype, переоформить весь класс через дефайны

        static inline const std::string tryToGetUninitializedField = "Try to get uninitialized field ";

        static DataBaseType stringToDataBaseType(const std::string& s);

        ArgsBase(int argc, char *argv[]);

        [[nodiscard]] DataBaseType dbType() const;

        [[nodiscard]] std::string dbName() const;

        ADD_ARGUMENT(int, test, 0)




        // TODO complete define

    protected:
        DataBaseType _dbType = DataBaseType::BadType;
        std::string _dbName;

        void parseDbType();

        void parseDbName();

        void checkDbType();

        void checkDbName();
    };
} // qtr



ABSL_FLAG(std::string, dbType, "",
          "Possible types: "
                  OPTION_NAME(QtrDrive) ", "
                  OPTION_NAME(QtrRam) ", "
                  OPTION_NAME(BingoNoSQL));

ABSL_FLAG(std::string, dbName, "",
          "Name of folders with data base's files");








