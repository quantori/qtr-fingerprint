#pragma once 

#include <gtest/gtest.h>

#include <filesystem>
#include <vector>

namespace qtr {

template<typename T>
void compareTwoVectors(const std::vector<T> &v1, const std::vector<T> &v2)
{
    EXPECT_EQ(v1.size(), v2.size());
    for(size_t i = 0; i < v1.size(); i++)
        EXPECT_EQ(v1.at(i), v2.at(i));
}

std::filesystem::path getDataDir();

} // namespace qtr