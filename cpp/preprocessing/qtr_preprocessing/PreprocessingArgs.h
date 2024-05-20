#pragma once

#include "ArgsBase.h"
#include <absl/flags/declare.h>
#include "PreprocessingType.h"

ABSL_DECLARE_FLAG(std::string, preprocessingType);

ABSL_DECLARE_FLAG(std::string, sourceDir);

ABSL_DECLARE_FLAG(std::string, destDir);

ABSL_DECLARE_FLAG(bool, properties);

ABSL_DECLARE_FLAG(std::string, molIdType);

ABSL_DECLARE_FLAG(std::string, fingerprintType);

namespace qtr {

    enum class MolIdType {
        BadType,
        SMILES,
        UID
    };

    enum class FingerprintType {
        BadType,
        Indigo,
        RDKit,
        Custom
    };

    class PreprocessingArgs : public ArgsBase {
        static inline const std::unordered_map<std::string, PreprocessingType> _strToPreprocessingType =
                {{FLAG_NAME(CSV), PreprocessingType::CSV},
                 {FLAG_NAME(SDF), PreprocessingType::SDF}};

        const static inline auto strToPreprocessingType = makeStringToEnumFunction(_strToPreprocessingType,
                                                                                   PreprocessingType::BadType);

        static inline const std::unordered_map<std::string, MolIdType> _strToMolIdType =
                {{FLAG_NAME(SMILES), MolIdType::SMILES},
                 {FLAG_NAME(UID),    MolIdType::UID}};

        const static inline auto strToMolIdType = makeStringToEnumFunction(_strToMolIdType, MolIdType::BadType);

        static inline const std::unordered_map<std::string, FingerprintType> _strToFingerprintType =
                {{FLAG_NAME(Indigo), FingerprintType::Indigo},
                 {FLAG_NAME(RDKit),  FingerprintType::RDKit},
                 {FLAG_NAME(Custom), FingerprintType::Custom}};

        const static inline auto strToFingerprintType = makeStringToEnumFunction(_strToFingerprintType,
                                                                                 FingerprintType::BadType);

    ADD_ARGUMENT_WITH_PARSER(PreprocessingType, preprocessingType, PreprocessingType::BadType, strToPreprocessingType);

    ADD_ARGUMENT(std::filesystem::path, sourceDir, "");

    ADD_ARGUMENT(std::filesystem::path, destDir, "");

    ADD_ARGUMENT(bool, properties, true);

    ADD_ARGUMENT_WITH_PARSER(FingerprintType, fingerprintType, FingerprintType::BadType, strToFingerprintType);


    ADD_ARGUMENT_WITH_PARSER(MolIdType, molIdType, MolIdType::BadType, strToMolIdType);

    public:
        PreprocessingArgs(int argc, char *argv[]) : ArgsBase(argc, argv) {
            parseAndCheck_preprocessingType();
            parseAndCheck_sourceDir();
            parseAndCheck_destDir();
            parseAndCheck_molIdType();
            if (preprocessingType() == PreprocessingType::SDF) {
            }
            if (preprocessingType() == PreprocessingType::CSV) {
                parse_properties();
                parseAndCheck_fingerprintType();
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

        [[nodiscard]] inline std::filesystem::path fingerprintLengthFile() const {
            return destDir() / "fingerprintLength";
        }
    };

} // qtr
