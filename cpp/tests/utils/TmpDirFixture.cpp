#include "TmpDirFixture.h"

const std::filesystem::path &TmpDirFixture::getTmpDir() const {
    return _dir;
}

void TmpDirFixture::SetUp() {
    _dir = qtr::DataPathManager::getTmpDataDir() / "tmpDirForTests";
    std::filesystem::remove_all(_dir);
    EXPECT_TRUE(std::filesystem::create_directory(_dir));
}

void TmpDirFixture::TearDown() {
    std::filesystem::remove_all(_dir);
}
