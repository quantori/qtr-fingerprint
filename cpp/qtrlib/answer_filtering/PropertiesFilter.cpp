#include "PropertiesFilter.h"
#include "glog/logging.h"

#include <algorithm>
#include <numeric>
#include <utility>
#include <cmath>

namespace qtr {
    PropertiesFilter::Bounds::Bounds() {
        for (size_t i = 0; i < PropertiesFilter::Properties::size(); i++) {
            minBounds[i] = -std::numeric_limits<property_t>::infinity();
            maxBounds[i] = std::numeric_limits<property_t>::infinity();
        }
    }

    bool PropertiesFilter::Bounds::Check(const PropertiesFilter::Properties &properties) const {
        for (size_t i = 0; i < PropertiesFilter::Properties::size(); i++) {
            if (!(minBounds[i] <= properties[i] && properties[i] <= maxBounds[i]) && !std::isnan(properties[i])) {
                return false;
            }
        }
        return true;
    }

    bool PropertiesFilter::operator()(CIDType id) {
        return _bounds->Check(_propertiesTable->operator[](id));
    }

    std::unique_ptr<AnswerFilter> PropertiesFilter::copy() {
        std::unique_ptr<PropertiesFilter> result = std::make_unique<PropertiesFilter>(_propertiesTable);
        result->setBounds(_bounds);
        return result;
    }

    PropertiesFilter::PropertiesFilter(std::shared_ptr<const std::vector<PropertiesFilter::Properties>> propertiesTable)
            : _propertiesTable(std::move(propertiesTable)), _bounds() {}

    PropertiesFilter::PropertiesFilter(std::shared_ptr<const std::vector<PropertiesFilter::Properties>> propertiesTable,
                                       PropertiesFilter::Bounds bounds) : PropertiesFilter(std::move(propertiesTable)) {
        setBounds(bounds);
    }

    void PropertiesFilter::setBounds(PropertiesFilter::Bounds bounds) {
        _bounds = std::make_shared<Bounds>(bounds);
    }

    void PropertiesFilter::setBounds(std::shared_ptr<const Bounds> bounds) {
        _bounds = std::move(bounds);
    }

    PropertiesFilter::property_t PropertiesFilter::Properties::operator[](size_t i) const {
        auto *arr = (property_t *) this;
        return arr[i];
    }

    PropertiesFilter::property_t &PropertiesFilter::Properties::operator[](size_t i) {
        auto *arr = (property_t *) this;
        return arr[i];
    }


} // qtr