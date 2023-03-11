#pragma once

#include <memory>
#include <filesystem>

#include "BallTreeTypes.h"

namespace qtr {
    class AnswerFilter {
    public:

        /// this class should be copied via method @copy
        AnswerFilter(const AnswerFilter &answerFilter) = delete;

        AnswerFilter() = default;

        virtual bool operator()(CIDType id) = 0;

        virtual std::unique_ptr<AnswerFilter> copy() = 0;

        virtual void initBallTreeLeaf(const std::filesystem::path &leafDirPath) {};

        virtual ~AnswerFilter() = default;
    };
}