#pragma once

#include <utility>
#include <vector>
#include <string>
#include <future>

#include "crow.h"

#include "modes/RunMode.h"
#include "QueryData.h"

namespace qtr {
    class WebMode : public RunMode {
    public:
        explicit WebMode(std::shared_ptr<SearchData> searchData);

        void run() override;

    private:
        crow::json::wvalue prepareResponse(QueryData<CIDType> &queryData, size_t minOffset, size_t maxOffset);
    };
} // namespace qtr
