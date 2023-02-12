#pragma once

#include "crow.h"

#include "Utils.h"
#include "io/BufferedReader.h"
#include "BallTreeSearchEngine.h"
#include "BallTreeRAMSearchEngine.h"
#include "fingerprint_table_io/FingerprintTableReader.h"
#include "IndigoSubstructureMatcher.h"
#include "IndigoQueryMolecule.h"
#include "string_table_io/StringTableReader.h"
#include "RunDbUtils.h"
#include "IdConverter.h"

#include "modes/RunMode.h"

namespace qtr {
    class WebMode : public RunMode {
    private:
        const qtr::BallTreeSearchEngine &_ballTree;
        std::shared_ptr<const SmilesTable> _smilesTable;
        const uint64_t _ansCount;
        const uint64_t _threadsCount;
        std::shared_ptr<const IdConverter> _idConverter;
        std::shared_ptr<const std::vector<PropertiesFilter::Properties>> _molPropertiesTable;

    public:
        WebMode(const qtr::BallTreeSearchEngine &ballTree, std::shared_ptr<const SmilesTable> smilesTable,
                qtr::TimeTicker &timeTicker, uint64_t ansCount, uint64_t threadsCount,
                std::shared_ptr<const IdConverter> idConverter,
                std::shared_ptr<const std::vector<PropertiesFilter::Properties>> molPropertiesTable);

        void run() override;

    private:
        crow::json::wvalue
        prepareResponse(BallTreeQueryData &queryData, size_t minOffset, size_t maxOffset);
    };
} // namespace qtr
