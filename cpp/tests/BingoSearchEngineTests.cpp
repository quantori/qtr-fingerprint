#include "common.h"
#include "SearchEngineFactory.h"

#include <gtest/gtest.h>

using namespace indigo_cpp;

TEST(BingoSearchEngineTest, Basic) 
{
    IndigoSessionPtr indigoSessionPtr = IndigoSession::create();
    SearchEnginePtr searchEnginePtr = SearchEngineFactory::create(indigoSessionPtr);

    testSearchEngine(searchEnginePtr, indigoSessionPtr);
}