#pragma once 

#include "SearchEngineInterface.h"

#include "IndigoSession.h"

#include <gtest/gtest.h>

#include <filesystem>
#include <vector>

void testSearchEngine(
    SearchEnginePtr searchEngine,
    indigo_cpp::IndigoSessionPtr indigoSession,
    const std::string &fileSdf,
    const std::string &fileQueries);

template<typename T>
void compareTwoVectors(const std::vector<T> &v1, const std::vector<T> &v2)
{
    EXPECT_EQ(v1.size(), v2.size());
    for(size_t i = 0; i < v1.size(); i++)
        EXPECT_EQ(v1.at(i), v2.at(i));
}

std::filesystem::path getDataDir();