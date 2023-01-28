#include "PropertiesFilter.h"

#include <algorithm>
#include <numeric>
#include <utility>

namespace qtr {
    bool PropertiesFilter::Bounds::Check(const PropertiesFilter::property_list_t &properties) const {
        for (size_t i = 0; i < propertiesCount; i++) {
            if (!(_minBounds[i] <= properties[i] && properties[i] <= _maxBounds[i]))
                return false;
        }
        return true;
    }

    PropertiesFilter::Bounds::Bounds() {
        std::fill(std::begin(_minBounds), std::end(_minBounds), std::numeric_limits<property_t>::min());
        std::fill(std::begin(_maxBounds), std::end(_maxBounds), std::numeric_limits<property_t>::max());
    }

    bool PropertiesFilter::operator()(CIDType id) {
        return _bounds->Check(_propertiesTable->operator[](id));
    }

    std::unique_ptr<AnswerFilter> PropertiesFilter::copy() {
        std::unique_ptr<PropertiesFilter> result = std::make_unique<PropertiesFilter>(_propertiesTable);
        result->setBounds(_bounds);
        return result;
    }

    PropertiesFilter::PropertiesFilter(std::shared_ptr<const std::vector<property_list_t>> propertiesTable)
            : _propertiesTable(std::move(propertiesTable)), _bounds() {}

    void PropertiesFilter::setBounds(PropertiesFilter::Bounds bounds) {
        _bounds = std::make_shared<Bounds>(bounds);
    }

    void PropertiesFilter::setBounds(std::shared_ptr<const Bounds> bounds) {
        _bounds = std::move(bounds);
    }

} // qtr