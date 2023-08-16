#pragma once

#include "BuildArgs.h"
#include "TimeMeasurer.h"

namespace qtr {

    class DatabaseBuilder {
    public:
        virtual void build(const BuildArgs &args, TimeMeasurer &timeMeasurer) = 0;

        virtual ~DatabaseBuilder() = default;
    };

} // qtr

