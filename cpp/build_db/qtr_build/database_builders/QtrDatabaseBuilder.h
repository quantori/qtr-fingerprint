#pragma once

#include "DatabaseBuilder.h"

namespace qtr {

    class QtrDatabaseBuilder : public DatabaseBuilder {
        void build(const qtr::BuildArgs &args) override;
    };

} // qtr
