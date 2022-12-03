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
#include "IdConverter.h"

#include "modes/RunMode.h"

namespace qtr {
    class WebMode : public RunMode {
    private:
        const qtr::BallTreeSearchEngine &_ballTree;
        const SmilesTable &_smilesTable;
        const uint64_t _ansCount;
        const uint64_t _startSearchDepth;
        IdConverter _idConverter;

    public:
        WebMode(const qtr::BallTreeSearchEngine &ballTree, const SmilesTable &smilesTable,
                qtr::TimeTicker &timeTicker, uint64_t ansCount, uint64_t startSearchDepth,
                std::filesystem::path &idToStringDirPath);

        void run() override;

    private:
        static crow::json::wvalue
        prepareResponse(BallTreeSearchEngine::QueryData &queryData, size_t minOffset, size_t maxOffset);
    };
} // namespace qtr
