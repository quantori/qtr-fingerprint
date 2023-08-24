#pragma once

#include "Preprocessor.h"

namespace qtr {
    class SDFPreprocessor : public Preprocessor {
        void run(const PreprocessingArgs &args) override;
    };
} // qtr

