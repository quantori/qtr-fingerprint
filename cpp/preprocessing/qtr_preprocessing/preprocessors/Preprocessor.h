#pragma once

#include "PreprocessingArgs.h"

namespace qtr {
    class Preprocessor {
    public:
        virtual void run(const PreprocessingArgs &args) = 0;

        virtual ~Preprocessor() = default;
    };
} // qtr

