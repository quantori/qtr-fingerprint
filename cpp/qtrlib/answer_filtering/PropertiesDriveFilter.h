#pragma once

#include <map>

#include "PropertiesFilter.h"

namespace qtr {

    class PropertiesDriveFilter : public PropertiesFilter {
    public:
        std::unique_ptr<AnswerFilter> copy() override;

        PropertiesDriveFilter() = default;

        explicit PropertiesDriveFilter(Bounds bounds);

    private:
        const Properties& getProperties(CIDType id) override;

        void initBallTreeLeaf(const std::filesystem::path& leafDirPath) override;

        std::map<CIDType, Properties> _propertiesTable;
    };

} // qtr
