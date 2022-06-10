#pragma once

#include "SearchEngineFactory.h"

#include "Common.h"

#include "IndigoSession.h"

#include <filesystem>

namespace qtr {

class SearchEngineTests {
public:
    SearchEngineTests(SearchEngineFactory::SearchEngineType searchEngineType)
        : _indigoSessionPtr(indigo_cpp::IndigoSession::create())
        , _searchEnginePtr(SearchEngineFactory::create(searchEngineType, _indigoSessionPtr))
        , _dataDir(getDataDir())
    { }

    void testPubchem10();
    void testPubchem100();
    void testPubchem994();

private:
    indigo_cpp::IndigoSessionPtr _indigoSessionPtr;
    SearchEnginePtr _searchEnginePtr;
    std::filesystem::path _dataDir;
};

} // namespace qtr