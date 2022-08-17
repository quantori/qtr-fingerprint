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

    std::filesystem::path DataPathManager::getTmpDataDir() {
        return _tmpDataDirPath;
    }

    void DataPathManager::init(const path &dataDirPath, const path &bigDataDirPath, const path& tmpDataDirPath)
    {
        if (dataDirPath.empty())
            throw std::invalid_argument("You have to set --data_dir_path flag");

        if (bigDataDirPath.empty())
            throw std::invalid_argument("You have to set --big_data_dir_path flag");

        if (tmpDataDirPath.empty())
            throw std::invalid_argument("You have to set --tmp_data_dir_path flag");

        _dataDirPath = dataDirPath;
        _bigDataDirPath = bigDataDirPath;
        _tmpDataDirPath = tmpDataDirPath;
    }
}