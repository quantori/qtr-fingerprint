#include "GTestEnvironment.h"
#include "DataPathManager.h"

namespace qtr
{
    GTestEnvironment::GTestEnvironment(const std::filesystem::path &bigDataDirPath) : _bigDataDirPath(bigDataDirPath)
    {
    }

    void GTestEnvironment::SetUp()
    {
        DataPathManager::init(_bigDataDirPath);
    }
}