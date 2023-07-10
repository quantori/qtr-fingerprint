#pragma once

#include <glog/logging.h>
#include <absl/flags/flag.h>
#include <absl/flags/parse.h>

#include <string>
#include <vector>
#include <filesystem>
#include <functional>
#include <map>

#include "Utils.h"

#define FLAG_NAME(x) #x

#define ADD_ARGUMENT_WITH_PARSER(type, name, emptyVal, parseFunction) \
    protected: \
        type _##name = emptyVal; \
        inline void check_##name() const { \
            if (checkEmpty(_##name, emptyVal)) { \
                LOG(INFO) << "Please specify " << #name << " option"; \
                exit(0); \
            } \
        } \
        inline void parse_##name() { \
            _##name = parseFunction(absl::GetFlag(FLAGS_##name)); \
        } \
        inline void parseAndCheck_##name() { \
            parse_##name(); \
            check_##name(); \
        }\
    public: \
        [[nodiscard]] inline type name() const { \
            if (checkEmpty(_##name, emptyVal)) { \
                LOG(ERROR) << "Try to get uninitialized  field " << #name; \
                exit(-1); \
            } \
            return _##name; \
        }


#define ADD_ARGUMENT(type, name, emptyVal) ADD_ARGUMENT_WITH_PARSER(type, name, emptyVal, std::identity())

namespace qtr {

    class ArgsBase {
    public:
        enum class DataBaseType {
            BadType,
            QtrDrive,
            QtrRam,
            BingoNoSQL
        };

        ArgsBase(int argc, char *argv[]);

    private:
        static inline const std::unordered_map<std::string, DataBaseType> _strToDataBaseType =
                {{"QtrDrive",   DataBaseType::QtrDrive},
                 {"QrtRam",     DataBaseType::QtrRam},
                 {"BingoNoSQL", DataBaseType::BingoNoSQL}};

    protected:
        static std::vector<std::filesystem::path> vecStrToVecPath(const std::vector<std::string> &v);

        static DataBaseType strToDataBaseType(const std::string &s);
    };
} // qtr








