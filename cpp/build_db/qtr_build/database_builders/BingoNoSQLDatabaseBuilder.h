#pragma once

#include "DatabaseBuilder.h"

namespace qtr {

    class BingoNoSQLDatabaseBuilder : public DatabaseBuilder {
    public:
        void build(const BuildArgs &args) override;

    };

} // qtr

