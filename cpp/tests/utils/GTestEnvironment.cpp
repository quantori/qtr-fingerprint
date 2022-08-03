#include "GTestEnvironment.h"

#include <utility>
#include "DataPathManager.h"

namespace qtr {
    GTestEnvironment::GTestEnvironment(std::filesystem::path dataDirPath,
                                       std::filesystem::path bigDataDirPath) : _dataDirPath(std::move(dataDirPath)),
                                                                               _bigDataDirPath(
                                                                                       std::move(bigDataDirPath)) {
    }

    void GTestEnvironment::SetUp() {
        DataPathManager::init(_dataDirPath, _bigDataDirPath);
    }
}