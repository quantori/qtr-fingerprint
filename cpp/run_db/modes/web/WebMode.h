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
        std::shared_ptr<const SearchData> _searchData;

    public:
        explicit WebMode(std::shared_ptr<const SearchData> searchData);

        void run() override;

    private:
        crow::json::wvalue
        prepareResponse(BallTreeQueryData &queryData, size_t minOffset, size_t maxOffset);
    };
} // namespace qtr
