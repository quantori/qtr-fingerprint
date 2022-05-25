#pragma once

#include <memory>

#include "SearchEngineInterface.h"
#include "BingoSearchEngine.h"

struct SearchEngineFactory {
    inline static SearchEnginePtr create(const indigo_cpp::IndigoSessionPtr &indigoSessionPtr) {
        return std::make_shared<BingoSearchEngine>(indigoSessionPtr);
    }
};