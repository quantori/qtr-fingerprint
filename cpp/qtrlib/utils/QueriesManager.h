#pragma once

#include <cstdint>
#include <vector>
#include <string>
#include <unordered_map>
#include <map>

#include "BallTreeSearchEngine.h"

namespace qtr {

    class QueriesManager {
    private:
        std::unordered_map<std::string, uint64_t> _smilesToQueryId;
        std::map<uint64_t, BallTreeSearchEngine::QueryData> _idToQuery;
        uint64_t _queryIdTicker = 0;
        BallTreeSearchEngine& _ballTree;


    public:

        struct QueryData {
            std::string querySmiles;
            size_t ansCount;
        };

        enum ResponseStatus {
            OK, // Full answer is returned
            PARTIAL, // Not all indexes exists
            CALCULATING, // Answers calculating still in progress
        };

        QueriesManager(BallTreeSearchEngine& ballTree);

        uint64_t createQuery(const std::string& query);

        void killQuery();

        std::vector<std::pair<ResponseStatus, std::vector<std::pair<uint64_t, std::string>>>>
        getAnswers(size_t beginIndex, size_t endIndex);
    };

} // qtr

