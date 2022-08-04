#include "GTestEnvironment.h"

#include <utility>
#include "DataPathManager.h"

namespace qtr {
    GTestEnvironment::GTestEnvironment(std::filesystem::path dataDirPath,
                                       std::filesystem::path bigDataDirPath,
                                       std::filesystem::path tmpDataDirPath) :
            _dataDirPath(std::move(dataDirPath)), _bigDataDirPath(std::move(bigDataDirPath)),
            _tmpDataDirPath(std::move(tmpDataDirPath)) {
    }

    void GTestEnvironment::SetUp() {
        DataPathManager::init(_dataDirPath, _bigDataDirPath, _tmpDataDirPath);
    }
}
