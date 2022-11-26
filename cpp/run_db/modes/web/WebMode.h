#pragma once

#include "crow.h"

#include "Utils.h"
#include "io/BufferedReader.h"
#include "BallTreeRAMSearchEngine.h"
#include "BallTreeNoChecksSearchEngine.h"
#include "fingerprint_table_io/FingerprintTableReader.h"
#include "IndigoSubstructureMatcher.h"
#include "IndigoQueryMolecule.h"
#include "smiles_table_io/SmilesTableReader.h"
#include "RunDbUtils.h"

#include "modes/RunMode.h"

namespace qtr {
    class WebMode : public RunMode {
    private:
        const qtr::BallTreeSearchEngine &ballTree;
        const SmilesTable &smilesTable;
        const uint64_t ansCount;
        const uint64_t startSearchDepth;
    public:
        WebMode(const qtr::BallTreeSearchEngine &ballTree, const SmilesTable &smilesTable,
                qtr::TimeTicker &timeTicker, uint64_t ansCount, uint64_t startSearchDepth);

        void run() override;

    private:
        crow::json::wvalue prepareResponse(const std::vector<uint64_t> &ids, size_t minOffset, size_t maxOffset);
    };
} // namespace qtr
