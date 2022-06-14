#include "Common.h"

#include "IndigoInChI.h"
#include "IndigoQueryMolecule.h"
#include "indigo.h"
#include "indigo-inchi.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

namespace qtr {

std::filesystem::path getDataDir()
{
    using namespace std::filesystem;
    const path currentDir = testing::UnitTest::GetInstance()->original_working_dir();
    const path dataDir = currentDir / path("../../data");
    return dataDir;
}

} // namespace qtr
