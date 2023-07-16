#pragma once

#include <memory>
#include "BallTreeTypes.h"

namespace qtr {

    template<typename T>
    class AnswerFilter {
    public:
        /// this class should be copied via method @copy
        AnswerFilter(const AnswerFilter &answerFilter) = delete;

        AnswerFilter() = default;

        virtual bool operator()(const T &value) = 0;

        virtual std::unique_ptr<AnswerFilter> copy() = 0;

        virtual ~AnswerFilter() = default;
    };

    using ByIdAnswerFilter = AnswerFilter<CIDType>;

} // qtr

