#include "PropertiesRamFilter.h"

using namespace std;

namespace qtr {
    unique_ptr <ByIdAnswerFilter> PropertiesRamFilter::copy() {
        unique_ptr<PropertiesRamFilter> result = make_unique<PropertiesRamFilter>(_propertiesTable);
        result->setBounds(_bounds);
        return result;
    }

    PropertiesRamFilter::PropertiesRamFilter(shared_ptr<const vector <Properties>> propertiesTable)
            : PropertiesFilter(), _propertiesTable(std::move(propertiesTable)) {}

    PropertiesRamFilter::PropertiesRamFilter(shared_ptr<const vector <Properties>> propertiesTable,
                                             PropertiesFilter::Bounds bounds) :
            PropertiesFilter(bounds), _propertiesTable(std::move(propertiesTable)) {}

    const PropertiesFilter::Properties &PropertiesRamFilter::getProperties(CIDType id) {
        return _propertiesTable->operator[](id);
    }

} // qtr