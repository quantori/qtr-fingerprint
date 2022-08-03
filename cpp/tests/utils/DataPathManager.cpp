#include "DataPathManager.h"

#include <gtest/gtest.h>

using namespace std::filesystem;

namespace qtr
{
    path DataPathManager::getDataDir()
    {
        const path currentDir = testing::UnitTest::GetInstance()->original_working_dir();
        const path dataDir = currentDir / path("../../../data");
        return dataDir;
    }

    path DataPathManager::getBigDataDir()
    {
        return _bigDataDirPath;
    }

    void DataPathManager::init(const path &bigDataDirPath)
    {
        _bigDataDirPath = bigDataDirPath;
    }
}