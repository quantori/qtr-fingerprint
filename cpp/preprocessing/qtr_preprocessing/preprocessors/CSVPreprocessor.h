#pragma once

#include "Preprocessor.h"

namespace qtr {
    class CSVPreprocessor : public Preprocessor {
        void run(const qtr::PreprocessingArgs &args) override;
    };
} // qtr

