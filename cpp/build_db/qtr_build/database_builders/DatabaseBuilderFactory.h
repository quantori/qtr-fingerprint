#pragma once

#include <memory>
#include "DatabaseBuilder.h"
#include "BuildArgs.h"

namespace qtr {

    class DatabaseBuilderFactory {
    public:
        static std::unique_ptr<DatabaseBuilder> create(DatabaseType databaseType);
    };

} // qtr
