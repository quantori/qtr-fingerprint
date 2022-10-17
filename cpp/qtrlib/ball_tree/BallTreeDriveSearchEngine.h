#pragma once

#include "BallTreeSearchEngine.h"

namespace qtr {

    class BallTreeDriveSearchEngine : public BallTreeSearchEngine {
    protected:
        virtual void searchInLeaf(size_t leafId, QueryData &queryData) const override;
    public:
        using BallTreeSearchEngine::BallTreeSearchEngine;
    };

} // qtr

