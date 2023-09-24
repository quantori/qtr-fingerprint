#include "PropertiesDriveFilter.h"
#include "properties_table_io/PropertiesTableReader.h"

using namespace std;

namespace qtr {
    PropertiesDriveFilter::PropertiesDriveFilter(PropertiesFilter::Bounds bounds) : PropertiesFilter(bounds) {}

    std::unique_ptr<ByIdAnswerFilter> PropertiesDriveFilter::copy() {
        unique_ptr<PropertiesDriveFilter> result = make_unique<PropertiesDriveFilter>();
        result->setBounds(_bounds);
        return result;
    }

    void PropertiesDriveFilter::initBallTreeLeaf(const filesystem::path &leafDirPath) {
        filesystem::path propertiesTablePath = leafDirPath / ("properties");
        PropertiesTableReader reader(propertiesTablePath);
        _propertiesTable.clear();
        _propertiesTable = map<CIDType, Properties>(reader.begin(), reader.end());
    }

    const PropertiesFilter::Properties &PropertiesDriveFilter::getProperties(CIDType id) {
        return _propertiesTable.at(id);
    }
} // qtr