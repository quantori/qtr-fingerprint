#pragma once

#include <vector>

#include "PropertiesFilter.h"

namespace qtr {

    class PropertiesRamFilter : public PropertiesFilter {
    public:
        std::unique_ptr<AnswerFilter> copy() override;

        explicit PropertiesRamFilter(std::shared_ptr<const std::vector<Properties>> propertiesTable);

        PropertiesRamFilter(std::shared_ptr<const std::vector<Properties>> propertiesTable, Bounds bounds);

    private:
        const Properties& getProperties(CIDType id) override;

        std::shared_ptr<const std::vector<Properties>> _propertiesTable;
    };

} // qtr
