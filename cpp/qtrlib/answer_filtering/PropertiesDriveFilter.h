#pragma once

#include <map>
#include <filesystem>

#include "PropertiesFilter.h"
#include "QtrBallTreeLeafInitMixin.h"

namespace qtr {

    class PropertiesDriveFilter : public PropertiesFilter, public QtrBallTreeLeafInitMixin {
    public:
        std::unique_ptr<AnswerFilter> copy() override;

        PropertiesDriveFilter() = default;

        explicit PropertiesDriveFilter(Bounds bounds);

        void initBallTreeLeaf(const std::filesystem::path& leafDirPath) override;

    private:
        const Properties& getProperties(CIDType id) override;

        std::map<CIDType, Properties> _propertiesTable;
    };

} // qtr
