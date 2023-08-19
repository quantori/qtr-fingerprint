#pragma once

#include <memory>

#include "Preprocessor.h"
#include "PreprocessingType.h"

namespace qtr {
    class PreprocessorFactory {
    public:
        static std::unique_ptr<Preprocessor> create(PreprocessingType type);
    };
} // qtr
