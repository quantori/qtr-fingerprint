#pragma once

#include "SearchEngineInterface.h"

#include <memory>

struct SearchEngineFactory {
    
    enum SearchEngineType {
        BINGO,
        EXHAUSTIVE
    };

    static SearchEnginePtr create(SearchEngineType searchEngineType, indigo_cpp::IndigoSessionPtr indigoSessionPtr);
};