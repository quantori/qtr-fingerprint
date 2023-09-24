#include "PreprocessorFactory.h"
#include "SDFPreprocessor.h"
#include "CSVPreprocessor.h"

namespace qtr {
    std::unique_ptr<Preprocessor> PreprocessorFactory::create(PreprocessingType type) {
        if (type == PreprocessingType::SDF) {
            return std::make_unique<SDFPreprocessor>();
        } else if (type == PreprocessingType::CSV) {
            return std::make_unique<CSVPreprocessor>();
        }
        throw std::logic_error("Preprocessor: Cannot create Preprocessor of given type");
    }
} // qtr