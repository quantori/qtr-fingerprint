#pragma once

#include "AnswerFilter.h"

namespace qtr {

    class CompoundFilter : public AnswerFilter {
    public:
        CompoundFilter(std::unique_ptr<AnswerFilter> filter1, std::unique_ptr<AnswerFilter> filter2);

        bool operator()(CIDType id) override;

        std::unique_ptr<AnswerFilter> copy() override;

    private:
        std::unique_ptr<AnswerFilter> _filter1;
        std::unique_ptr<AnswerFilter> _filter2;
    };

} // qtr
