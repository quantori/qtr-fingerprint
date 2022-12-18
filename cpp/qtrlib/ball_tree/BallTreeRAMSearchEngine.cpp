#include "BallTreeRAMSearchEngine.h"

#include "fingerprint_table_io/FingerprintTableReader.h"

#include "IndigoQueryMolecule.h"
#include "IndigoSubstructureMatcher.h"

namespace qtr {

    std::map<std::string, std::vector<std::pair<size_t, std::filesystem::path>>>
    BallTreeRAMSearchEngine::groupedByDriveLeafFiles() const {
        std::map<std::string, std::vector<std::pair<size_t, std::filesystem::path>>> result;
        for (size_t leafId: getLeafIds()) {
            auto leafFilePath = getLeafFile(leafId);
            result[leafFilePath.parent_path().parent_path().string()].emplace_back(leafId, leafFilePath);
        }
        return result;
    }

    void
    BallTreeRAMSearchEngine::loadLeafFiles(const std::vector<std::pair<size_t, std::filesystem::path>> &leafsList) {
        LOG(INFO) << "Start loading leaf content from drive " << leafsList.front().second.parent_path().parent_path();
        for (auto &[leafId, leafFilePath]: leafsList) {
            FingerprintTableReader reader(leafFilePath);
            std::copy(reader.begin(), reader.end(), std::back_inserter(_buckets[leafNumberById(leafId)]));
        }
        LOG(INFO) << "Finish loading leaf content from drive " << leafsList.front().second.parent_path().parent_path();
    }

    void BallTreeRAMSearchEngine::searchInLeafIds(const indigo_cpp::IndigoSessionPtr &indigoSessionPtr,
                                                  const indigo_cpp::IndigoQueryMolecule &queryMol,
                                                  size_t leafId, BallTreeSearchEngine::QueryData &queryData) const {
        for (const auto &[id, fingerprint]: _buckets[leafNumberById((leafId))]) {
            if (queryData.query <= fingerprint)
            {
                if (queryData.isTerminate)
                    break;
                const auto &ansSmiles = queryData.smilesTable.at(id);
                try {
                    auto candidateMol = indigoSessionPtr->loadMolecule(ansSmiles);
                    candidateMol.aromatize();
                    auto matcher = indigoSessionPtr->substructureMatcher(candidateMol);
                    if (bool(indigoMatch(matcher.id(), queryMol.id())))
                        queryData.addAnswer(id);
                }
                catch (std::exception &e) {
                    LOG(ERROR) << "Error while filtering answer. "
                                  "Query: " << queryData.querySmiles << ", candidate: " << id << " " << ansSmiles << ", error: "
                               << e.what();
                }
            }
        }
    }

    void BallTreeRAMSearchEngine::searchInLeaf(size_t leafId, BallTreeSearchEngine::QueryData &queryData) const {
        for (const auto &[id, fingerprint]: _buckets[leafNumberById((leafId))]) {
            if (queryData.query <= fingerprint)
                putAnswer(id, queryData);
        }
    }

} // qtr
