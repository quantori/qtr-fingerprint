#include "RunDbUtils.h"

#include "smiles_table_io/SmilesTableReader.h"
#include "BallTreeQueryData.h"
#include "IndigoQueryMolecule.h"
#include "IndigoSubstructureMatcher.h"
#include "IndigoFilter.h"
#include "PropertiesFilter.h"
#include "CompoundFilter.h"

using namespace std;

namespace qtr {

    std::pair<bool, std::unique_ptr<BallTreeQueryData>>
    doSearch(const std::string &querySmiles, const qtr::BallTreeSearchEngine &ballTree,
             const std::shared_ptr<const SmilesTable> &smilesTable, size_t ansCount, size_t threadsCount,
             const PropertiesFilter::Bounds &queryBounds,
             const shared_ptr<const vector<PropertiesFilter::Properties>> &propertiesTable) {
        qtr::IndigoFingerprint fingerprint;
        try {
            fingerprint = qtr::indigoFingerprintFromSmiles(querySmiles);
        }
        catch (exception &exception) {
            LOG(WARNING) << "Skip query:" << exception.what();
            return {true, unique_ptr<BallTreeQueryData>(nullptr)};
        }
        LOG(INFO) << "Start search: " << querySmiles;
        auto querySmilesPtr = make_shared<string>(querySmiles);
        auto indigoFilter = std::make_unique<IndigoFilter>(smilesTable, querySmilesPtr);
        auto propertiesFilter = std::make_unique<PropertiesFilter>(propertiesTable, queryBounds);
        auto filter = std::make_unique<CompoundFilter>(std::move(propertiesFilter), std::move(indigoFilter));

        auto queryData = make_unique<BallTreeQueryData>(ansCount, fingerprint, std::move(filter));
        ballTree.search(*queryData, threadsCount);
        return {false, std::move(queryData)};
    }

    std::shared_ptr<SmilesTable>
    loadSmilesTable(const filesystem::path &smilesTablePath, const HuffmanCoder &huffmanCoder) {
        LOG(INFO) << "Start smiles table loading";
        HuffmanSmilesTable::Builder builder(huffmanCoder);
        for (const auto &pair: SmilesTableReader(smilesTablePath)) {
            builder += pair;
        }
        LOG(INFO) << "Finish smiles table loading";
        return builder.buildPtr();
    }

} // namespace qtr

