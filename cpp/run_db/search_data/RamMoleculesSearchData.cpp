#include "RamMoleculesSearchData.h"

namespace qtr {
    DbType RamMoleculesSearchData::getDbType() const {
        return DbType::InRamMolecules;
    }

    RamMoleculesSearchData::RamMoleculesSearchData(std::shared_ptr<const BallTreeSearchEngine> ballTree,
                                                   std::shared_ptr<const IdConverter> idConverter,
                                                   TimeTicker &timeTicker, size_t ansCount, size_t threadCount,
                                                   std::shared_ptr<const std::map<CIDType, indigo_cpp::IndigoMolecule>> moleculesTable,
                                                   std::shared_ptr<const std::vector<PropertiesFilter::Properties>> propertiesTable)
            :
            SearchData(std::move(ballTree), std::move(idConverter), timeTicker, ansCount, threadCount),
            moleculesTable(std::move(moleculesTable)), propertiesTable(std::move(propertiesTable)) {}
} // qtr