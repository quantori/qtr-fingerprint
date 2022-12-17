#pragma once

#include "AnswerFilter.h"

namespace qtr {

    class AlwaysTrueFilter : public AnswerFilter {
    public:
        AlwaysTrueFilter() = default;

        inline bool operator()(CIDType id) override {
            return true;
        }

        inline std::unique_ptr<AnswerFilter> copy() override {
            return std::make_unique<AlwaysTrueFilter>();
        }
    };

} // qtr
