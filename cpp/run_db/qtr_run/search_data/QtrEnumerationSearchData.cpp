#include "QtrEnumerationSearchData.h"

#include <utility>

using namespace std;

namespace qtr {
    unique_ptr <QueryData<CIDType>>
    QtrEnumerationSearchData::search(const SearchData::Query &query, const PropertiesFilter::Bounds &) {
        Fingerprint fingerprint = query.getFingerprint();
        auto queryData = make_unique<QueryDataWithFingerprint>(ansCount, timeLimit, fingerprint,
                                                               make_unique<AlwaysTrueFilter<CIDType>>());

        auto searchFunction = [this](QueryDataWithFingerprint& queryData) {
            auto& queryFingerprint = queryData.getQueryFingerprint();
            assert(queryFingerprint.size() == 1024);
            std::vector<CIDType> answers;
            for (size_t i = 0; i < _allFingerprints->size() && i < ansCount; i++) {
                auto& fingerprint = _allFingerprints->at(i);
                if (queryFingerprint <= fingerprint)
                    answers.emplace_back(i);
                if (i & ((1 << 16) - 1)) { // magic constant
                    queryData.addAnswers(answers);
                    answers.clear();
                    if (queryData.checkTimeOut())
                        break;
                }
            }
            queryData.addAnswers(answers);
            if (queryData.checkTimeOut()) {
                LOG(INFO) << "Search stopped due to timeout";
            }
            queryData.tagFinishTask();
        };

        auto task = async(std::launch::async, searchFunction, std::ref(*queryData));
        queryData->addTask(std::move(task));
        return std::move(queryData);
    }

    QtrEnumerationSearchData::QtrEnumerationSearchData(std::shared_ptr<std::vector<Fingerprint>> allFingerprints,
                                                       size_t ansCount, size_t threadCount, double timeLimit,
                                                       bool verificationStage) : SearchData(ansCount, threadCount, timeLimit, verificationStage),
                                                                                 _allFingerprints(std::move(allFingerprints)) {
        assert(threadCount == 1);
        assert(!verificationStage);
    }
} // qtr