#pragma once

#include "AnswerFilter.h"

namespace qtr {

    template<typename T>
    class AlwaysTrueFilter : public AnswerFilter<T> {
    public:
        AlwaysTrueFilter() = default;

        inline bool operator()(const T &) override {
            return true;
        }

        inline std::unique_ptr<AnswerFilter<T>> copy() override {
            return std::make_unique<AlwaysTrueFilter>();
        }
    };

} // qtr
