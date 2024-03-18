#include "QtrEnumerationSearchData.h"

#include <utility>

using namespace std;

namespace qtr {
    unique_ptr <QueryData<CIDType>>
    QtrEnumerationSearchData::search(const SearchData::Query &query, const PropertiesFilter::Bounds &) {
        Fingerprint fingerprint = query.getFingerprint();
        auto queryData = make_unique<QueryDataWithFingerprint>(ansCount, timeLimit, fingerprint,
                                                               make_unique<AlwaysTrueFilter<CIDType>>());

        assert(!verificationStage);

        auto searchFunction = [this](QueryDataWithFingerprint &queryData) {
            auto &queryFingerprint = queryData.getQueryFingerprint();
            for (auto &leafId: ballTree->getLeafIds()) {
                if (queryData.checkTimeOut()) {
                    LOG(INFO) << "Search stopped due to timeout";
                    break;
                }
                if (queryData.getCurrentAnswersCount() >= ansCount) {
                    break;
                }
                std::vector<CIDType> answers;
                for (auto &[id, fp]: ballTree->getLeafContent(leafId)) {
                    if (queryFingerprint <= fp)
                        answers.emplace_back(id);
                }
                queryData.addAnswers(answers);
            }
            queryData.tagFinishTask();
        };

        auto task = async(std::launch::async, searchFunction, std::ref(*queryData));
        queryData->addTask(std::move(task));
        return std::move(queryData);
    }

    QtrEnumerationSearchData::QtrEnumerationSearchData(QtrRamSearchData searchData) : QtrRamSearchData(
            std::move(searchData)) {
    }
} // qtr