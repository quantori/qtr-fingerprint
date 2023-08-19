#pragma once

#include "ArgsBase.h"
#include <absl/flags/declare.h>
#include "PreprocessingType.h"

ABSL_DECLARE_FLAG(std::string, preprocessingType);

ABSL_DECLARE_FLAG(std::string, sourceDir);

ABSL_DECLARE_FLAG(std::string, destDir);

ABSL_DECLARE_FLAG(std::string, destFilesType);


namespace qtr {

    enum class DestinationType {
        BadType,
        RawBucket,
        Tables
    };

    class PreprocessingArgs : public ArgsBase {
        static inline const std::unordered_map<std::string, DestinationType> _strToDestinationType =
                {{"RawBucket", DestinationType::RawBucket},
                 {"Tables",    DestinationType::Tables}};

        const static inline auto strToDestinationType = makeStringToEnumFunction(_strToDestinationType,
                                                                                 DestinationType::BadType);

        static inline const std::unordered_map<std::string, PreprocessingType> _strToPreprocessingType =
                {{"csv", PreprocessingType::CSV},
                 {"sdf", PreprocessingType::SDF}};

        const static inline auto strToPreprocessingType = makeStringToEnumFunction(_strToPreprocessingType,
                                                                                   PreprocessingType::BadType);

    ADD_ARGUMENT_WITH_PARSER(PreprocessingType, preprocessingType, PreprocessingType::BadType, strToPreprocessingType);

    ADD_ARGUMENT(std::filesystem::path, sourceDir, "");

    ADD_ARGUMENT(std::filesystem::path, destDir, "");

    ADD_ARGUMENT_WITH_PARSER(DestinationType, destFilesType, DestinationType::BadType, strToDestinationType);

    public:
        PreprocessingArgs(int argc, char *argv[]) : ArgsBase(argc, argv) {
            parseAndCheck_preprocessingType();
            parseAndCheck_sourceDir();
            parseAndCheck_destDir();
            if (preprocessingType() == PreprocessingType::SDF) {
                parse_destFilesType();
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
