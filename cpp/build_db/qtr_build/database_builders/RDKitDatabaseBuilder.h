#pragma once

#include "DatabaseBuilder.h"

namespace qtr {


    class RDKitDatabaseBuilder : public DatabaseBuilder {
    public:
        void build(const BuildArgs &args) override;
    };

}
