#pragma once

#include <filesystem>

namespace qtr {

    class QtrBallTreeLeafInitMixin {
    public:
        virtual void initBallTreeLeaf(const std::filesystem::path &leafDir) = 0;
    };

    template<typename AnsType>
    void tryInitBallTreeLeaf(AnswerFilter<AnsType> &filter, const std::filesystem::path &leafDir) {
        auto leafInit = dynamic_cast<QtrBallTreeLeafInitMixin *> (&filter);
        if (leafInit != nullptr)
            leafInit->initBallTreeLeaf(leafDir);
    }

} // qtr

