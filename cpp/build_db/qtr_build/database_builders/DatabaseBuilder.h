#pragma once

#include "BuildArgs.h"

namespace qtr {

    class DatabaseBuilder {
    public:
        virtual void build(const BuildArgs &args) = 0;

        virtual ~DatabaseBuilder() = default;
    };

} // qtr

