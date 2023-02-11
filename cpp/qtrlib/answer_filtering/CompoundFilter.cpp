#include "CompoundFilter.h"

namespace qtr {
    bool CompoundFilter::operator()(CIDType id) {
        return _filter1->operator()(id) && _filter2->operator()(id);
    }

    std::unique_ptr<AnswerFilter> CompoundFilter::copy() {
        return std::make_unique<CompoundFilter>(_filter1->copy(), _filter2->copy());
    }

    CompoundFilter::CompoundFilter(std::unique_ptr<AnswerFilter> filter1, std::unique_ptr<AnswerFilter> filter2)
            : _filter1(std::move(filter1)), _filter2(std::move(filter2)) {}
} // qtr