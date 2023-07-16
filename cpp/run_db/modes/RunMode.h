#pragma once

#include <utility>

#include "search_data/SearchData.h"

namespace qtr {
    class RunMode {
    protected:
        std::shared_ptr<SearchData> _searchData;

    public:
        explicit inline RunMode(std::shared_ptr<SearchData> searchData) : _searchData(std::move(searchData)) {}

        virtual void run() = 0;

        virtual ~RunMode() = default;
    };
}