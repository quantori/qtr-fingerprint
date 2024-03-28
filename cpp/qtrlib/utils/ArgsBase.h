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
#include "DatabaseType.h"

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
            _##name = parseFunction(absl::GetFlag(FLAGS_##name));     \
            LOG(INFO) << #name << ": " << toString(absl::GetFlag(FLAGS_##name)); \
        } \
        inline void parseAndCheck_##name() { \
            parse_##name(); \
            check_##name(); \
        }\
    public: \
        [[nodiscard]] inline type name() const { \
            return _##name; \
        } \
        [[nodiscard]] inline bool isProvided_##name() const { \
            return checkEmpty(_##name, emptyVal); \
        } \


#define ADD_ARGUMENT(type, name, emptyVal) ADD_ARGUMENT_WITH_PARSER(type, name, emptyVal, std::identity())

namespace qtr {

    class ArgsBase {
    public:
        ArgsBase(int argc, char *argv[]);

    private:
        static inline const std::unordered_map<std::string, DatabaseType> _strToDataBaseType =
                {{"QtrDrive",    DatabaseType::QtrDrive},
                 {"QtrRam",      DatabaseType::QtrRam},
                 {"BingoNoSQL",  DatabaseType::BingoNoSQL},
                 {"QtrEnumeration", DatabaseType::QtrEnumeration},};

    protected:
        static std::vector<std::filesystem::path> vecStrToVecPath(const std::vector<std::string> &v);

        const static inline auto strToDataBaseType = makeStringToEnumFunction(_strToDataBaseType,
                                                                              DatabaseType::BadType);

    };
} // qtr








