#include "RunDbUtils.h"

#include "smiles_table_io/SmilesTableReader.h"
#include "BallTreeQueryData.h"
#include "IndigoQueryMolecule.h"
#include "IndigoSubstructureMatcher.h"
#include "answer_filtering/IndigoFilter.h"

using namespace std;

namespace qtr {

    pair<bool, std::unique_ptr<BallTreeQueryData>>
    doSearch(const string &querySmiles, const qtr::BallTreeSearchEngine &ballTree,
             const shared_ptr<const SmilesTable> &smilesTable, size_t ansCount, size_t threadsCount) {
        qtr::IndigoFingerprint fingerprint;
        try {
            fingerprint = qtr::indigoFingerprintFromSmiles(querySmiles);
        }
        catch (exception &exception) {
            cout << "skip query:" << exception.what() << endl;
            return {true, unique_ptr<BallTreeQueryData>(nullptr)};
        }
        auto querySmilesPtr = make_shared<string>(querySmiles);
        auto queryData = make_unique<BallTreeQueryData>(ansCount, fingerprint,
                                                        std::make_unique<IndigoFilter>(smilesTable, querySmilesPtr));
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

