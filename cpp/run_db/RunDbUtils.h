#pragma once

using SmilesTable = std::unordered_map<uint64_t, std::string>;

namespace qtr {
    std::pair<bool, std::vector<uint64_t>>
    inline doSearch(const std::string &querySmiles, const qtr::BallTreeSearchEngine &ballTree,
                    const SmilesTable &smilesTable, uint64_t ansCount, uint64_t startSearchDepth) {
        qtr::IndigoFingerprint fingerprint;
        try {
            fingerprint = qtr::IndigoFingerprintFromSmiles(querySmiles);
        }
        catch (std::exception &exception) {
            std::cout << "skip query:" << exception.what() << std::endl;
            return {true, {}};
        }

        auto filter = [&smilesTable, &querySmiles](size_t ansId) {
            const auto &ansSmiles = smilesTable.at(ansId);
            auto indigoSessionPtr = indigo_cpp::IndigoSession::create();
            auto queryMol = indigoSessionPtr->loadQueryMolecule(querySmiles);
            queryMol.aromatize();
            try {
                auto candidateMol = indigoSessionPtr->loadMolecule(ansSmiles);
                candidateMol.aromatize();
                auto matcher = indigoSessionPtr->substructureMatcher(candidateMol);
                return bool(indigoMatch(matcher.id(), queryMol.id()));
            }
            catch (std::exception &e) {
                LOG(ERROR) << "Error while filtering answer. "
                              "Query: " << querySmiles << ", candidate: " << ansId << " " << ansSmiles << ", error: "
                           << e.what();
                return false;
            }
        };
        auto candidateIndexes = ballTree.search(fingerprint, ansCount, startSearchDepth, filter);
        LOG(INFO) << "found answers: " << candidateIndexes.size();
        return {false, candidateIndexes};
    }

    inline void
    loadSmilesTable(SmilesTable &smilesTable, const std::filesystem::path &smilesTablePath) {
        LOG(INFO) << "Start smiles table loading";
        for (const auto &value: qtr::SmilesTableReader(smilesTablePath)) {
            smilesTable[value.first] = value.second;
        }
        LOG(INFO) << "Finish smiles table loading";
    }
} // namespace qtr