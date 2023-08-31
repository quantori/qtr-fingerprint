#pragma once

#include "ArgsBase.h"
#include <absl/flags/declare.h>
#include "PreprocessingType.h"

ABSL_DECLARE_FLAG(std::string, preprocessingType);

ABSL_DECLARE_FLAG(std::string, sourceDir);

ABSL_DECLARE_FLAG(std::string, destDir);

ABSL_DECLARE_FLAG(std::string, targetFilesType);

ABSL_DECLARE_FLAG(bool, properties);

namespace qtr {

    enum class TargetType {
        BadType,
        RawBucket,
        Tables
    };

    class PreprocessingArgs : public ArgsBase {
        static inline const std::unordered_map<std::string, TargetType> _strToDestinationType =
                {{FLAG_NAME(RawBucket), TargetType::RawBucket},
                 {FLAG_NAME(Tables),    TargetType::Tables}};

        const static inline auto strToTargetType = makeStringToEnumFunction(_strToDestinationType,
                                                                            TargetType::BadType);

        static inline const std::unordered_map<std::string, PreprocessingType> _strToPreprocessingType =
                {{FLAG_NAME(CSV), PreprocessingType::CSV},
                 {FLAG_NAME(SDF), PreprocessingType::SDF}};

        const static inline auto strToPreprocessingType = makeStringToEnumFunction(_strToPreprocessingType,
                                                                                   PreprocessingType::BadType);

    ADD_ARGUMENT_WITH_PARSER(PreprocessingType, preprocessingType, PreprocessingType::BadType, strToPreprocessingType);

    ADD_ARGUMENT(std::filesystem::path, sourceDir, "");

    ADD_ARGUMENT(std::filesystem::path, destDir, "");

    ADD_ARGUMENT_WITH_PARSER(TargetType, targetFilesType, TargetType::BadType, strToTargetType);

    ADD_ARGUMENT(bool, properties, true);

    public:
        PreprocessingArgs(int argc, char *argv[]) : ArgsBase(argc, argv) {
            parseAndCheck_preprocessingType();
            parseAndCheck_sourceDir();
            parseAndCheck_destDir();
            if (preprocessingType() == PreprocessingType::SDF) {
                parse_targetFilesType();
            }
            if (preprocessingType() == PreprocessingType::CSV) {
                parse_properties();
            }
        }

        [[nodiscard]] inline std::filesystem::path smilesTables() const {
            return destDir() / "smilesTables";
        }

        [[nodiscard]] inline std::filesystem::path fingerprintTables() const {
            return destDir() / "fingerprintTables";
        }

        [[nodiscard]] inline std::filesystem::path idToStringTables() const {
            return destDir() / "idToStringTables";
        }

        [[nodiscard]] inline std::filesystem::path propertyTables() const {
            return destDir() / "propertyTables";
        }
    };

} // qtr
