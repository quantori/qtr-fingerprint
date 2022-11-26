#pragma once

#include "BallTreeSearchEngine.h"

namespace qtr {

    class BallTreeNoChecksSearchEngine : public BallTreeSearchEngine {
    private:
        std::vector<std::vector<uint64_t>> _leafData;

        void loadLeafData();

        void searchInLeaf(size_t leafId, QueryData& queryData) const override;

    public:
        template<typename BinaryReader>
        BallTreeNoChecksSearchEngine(BinaryReader &nodesReader,
                                     std::vector<std::filesystem::path> dataDirectories);
    };

    template<typename BinaryReader>
    BallTreeNoChecksSearchEngine::BallTreeNoChecksSearchEngine(BinaryReader &nodesReader,
                                                               std::vector<std::filesystem::path> dataDirectories)
            :BallTreeSearchEngine(nodesReader, dataDirectories) {
        loadLeafData();
    }

} // qtr
