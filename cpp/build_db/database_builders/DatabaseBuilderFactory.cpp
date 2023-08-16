#include "DatabaseBuilderFactory.h"

#include "QtrDriveDatabaseBuilder.h"
#include "QtrRamDatabaseBuilder.h"
#include "BingoNoSQLDatabaseBuilder.h"

using namespace std;

namespace qtr {
    unique_ptr<DatabaseBuilder> DatabaseBuilderFactory::create(DatabaseType databaseType) {
        if (databaseType == DatabaseType::QtrDrive)
            return make_unique<QtrDriveDatabaseBuilder>();
        else if (databaseType == DatabaseType::QtrRam)
            return make_unique<QtrRamDatabaseBuilder>();
        else if (databaseType == DatabaseType::BingoNoSQL)
            return make_unique<BingoNoSQLDatabaseBuilder>();

        throw logic_error("DatabaseBuilderFactory: Cannot create DatabaseBuilder of given type");
    }
} // qtr