#include "DataPathManager.h"

#include <gtest/gtest.h>

using namespace std::filesystem;

namespace qtr
{
    path DataPathManager::getDataDir()
    {
        return _dataDirPath;
    }

    path DataPathManager::getBigDataDir()
    {
        return _bigDataDirPath;
    }

    void DataPathManager::init(const path &dataDirPath, const path &bigDataDirPath)
    {
        _dataDirPath = dataDirPath;
        _bigDataDirPath = bigDataDirPath;
    }
}