#pragma once

#include <map>

#include "SearchData.h"


namespace qtr {

    class RamMoleculesSearchData : public SearchData {
    public:
        RamMoleculesSearchData(std::shared_ptr<const BallTreeSearchEngine> ballTree,
        std::shared_ptr<const IdConverter> idConverter, TimeTicker &timeTicker, size_t ansCount,
                size_t threadCount, std::shared_ptr<const std::map<CIDType, indigo_cpp::IndigoMolecule>> moleculesTable,
                std::shared_ptr<const std::vector<PropertiesFilter::Properties>> propertiesTable);

        [[nodiscard]] DbType getDbType() const override;

        std::shared_ptr<const std::map<CIDType, indigo_cpp::IndigoMolecule>> moleculesTable;
        std::shared_ptr<const std::vector<PropertiesFilter::Properties>> propertiesTable;
    };

} // qtr

