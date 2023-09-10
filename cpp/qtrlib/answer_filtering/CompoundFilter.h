#pragma once

#include "AnswerFilter.h"
#include "QtrBallTreeLeafInitMixin.h"

namespace qtr {

    template<typename T>
    class CompoundFilter : public AnswerFilter<T>, public QtrBallTreeLeafInitMixin {
    public:
        CompoundFilter(std::unique_ptr<AnswerFilter<T>> filter1, std::unique_ptr<AnswerFilter<T>> filter2);

        bool operator()(const T &value) override;

        std::unique_ptr<AnswerFilter<T>> copy() override;

        void initBallTreeLeaf(const std::filesystem::path &leafDir) override;

    private:
        std::unique_ptr<AnswerFilter<T>> _filter1;
        std::unique_ptr<AnswerFilter<T>> _filter2;
    };

    template<typename T>
    void CompoundFilter<T>::initBallTreeLeaf(const std::filesystem::path &leafDir) {
        tryInitBallTreeLeaf(*_filter1, leafDir);
        tryInitBallTreeLeaf(*_filter2, leafDir);
    }

    template<typename T>
    bool CompoundFilter<T>::operator()(const T &value) {
        ProfileScope("Compound filter");
        return _filter1->operator()(value) && _filter2->operator()(value);
    }

    template<typename T>
    std::unique_ptr<AnswerFilter<T>> CompoundFilter<T>::copy() {
        return std::make_unique<CompoundFilter>(_filter1->copy(), _filter2->copy());
    }

    template<typename T>
    CompoundFilter<T>::CompoundFilter(std::unique_ptr<AnswerFilter<T>> filter1,
                                      std::unique_ptr<AnswerFilter<T>> filter2)
            : _filter1(std::move(filter1)), _filter2(std::move(filter2)) {}

} // qtr
