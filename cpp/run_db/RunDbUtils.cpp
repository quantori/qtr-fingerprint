#include "RunDbUtils.h"

#include "smiles_table_io/SmilesTableReader.h"

#include "IndigoQueryMolecule.h"
#include "IndigoSubstructureMatcher.h"

using namespace std;

namespace qtr {

    pair<bool, vector<uint64_t>>
    doSearch(const string &querySmiles, const qtr::BallTreeSearchEngine &ballTree,
             const SmilesTable &smilesTable, uint64_t ansCount, uint64_t startSearchDepth) {
        qtr::IndigoFingerprint fingerprint;
        try {
            fingerprint = qtr::IndigoFingerprintFromSmiles(querySmiles);
        }
        catch (exception &exception) {
            cout << "skip query:" << exception.what() << endl;
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
            catch (exception &e) {
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

    HuffmanCoder buildHuffmanCoder(const filesystem::path &smilesTablePath) {
        LOG(INFO) << "Start huffman coder building";
        qtr::HuffmanCoder::Builder huffmanBuilder;
        for (const auto &[_, smiles]: qtr::SmilesTableReader(smilesTablePath)) {
            huffmanBuilder += smiles;
        }
        LOG(INFO) << "Finish huffman coder building";
        return huffmanBuilder.build();
    }

    SmilesTable loadSmilesTable(const filesystem::path &smilesTablePath, const HuffmanCoder &huffmanCoder) {
        LOG(INFO) << "Start smiles table loading";
        SmilesTable::Builder builder(huffmanCoder);
        for (const auto &pair: qtr::SmilesTableReader(smilesTablePath)) {
            builder += pair;
        }
        LOG(INFO) << "Finish smiles table loading";
        return builder.build();
    }

    pair<HuffmanCoder, SmilesTable> loadCoderAndTable(const filesystem::path &smilesTablePath) {
        auto coder = buildHuffmanCoder(smilesTablePath);
        auto table = loadSmilesTable(smilesTablePath, coder);
        return {coder, table};
    }

} // namespace qtr

