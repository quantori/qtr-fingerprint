#include "PropertiesFilter.h"
#include "glog/logging.h"

#include <algorithm>
#include <numeric>
#include <utility>
#include <cmath>

namespace qtr {
    PropertiesFilter::Bounds::Bounds() {
        for (size_t i = 0; i < minBounds.size(); i++) {
            minBounds[i] = -std::numeric_limits<property_t>::infinity();
            maxBounds[i] = std::numeric_limits<property_t>::infinity();
        }
    }

    bool PropertiesFilter::Bounds::Check(const PropertiesFilter::Properties &properties) const {
        const static auto eps = (float) 1e-6;
        for (size_t i = 0; i < properties.size(); i++) {
            if (!std::isnan(properties[i]) &&
                (properties[i] < minBounds[i] - eps || properties[i] > maxBounds[i] + eps)) {
                return false;
            }
        }
        return true;
    }

    bool PropertiesFilter::operator()(CIDType id) {
        return _bounds->Check(getProperties(id));
    }

    void PropertiesFilter::setBounds(PropertiesFilter::Bounds bounds) {
        _bounds = std::make_shared<Bounds>(bounds);
    }

    void PropertiesFilter::setBounds(std::shared_ptr<const Bounds> bounds) {
        _bounds = std::move(bounds);
    }

    PropertiesFilter::PropertiesFilter() : _bounds() {}

    PropertiesFilter::PropertiesFilter(PropertiesFilter::Bounds bounds) : PropertiesFilter() {
        setBounds(bounds);
    }

    PropertiesFilter::property_t PropertiesFilter::Properties::operator[](size_t i) const {
        return _propertiesArr[i];
    }

    PropertiesFilter::property_t &PropertiesFilter::Properties::operator[](size_t i) {
        return _propertiesArr[i];
    }

    PropertiesFilter::property_t PropertiesFilter::Properties::pubchemComponentCount() const {
        return _propertiesArr[ComponentCount];
    }

    PropertiesFilter::property_t PropertiesFilter::Properties::pubchemXlogp3() const {
        return _propertiesArr[Xlogp3];
    }

    PropertiesFilter::property_t PropertiesFilter::Properties::pubchemAtomUdefStereoCount() const {
        return _propertiesArr[AtomUdefStereoCount];
    }

    PropertiesFilter::property_t PropertiesFilter::Properties::pubchemHeavyAtomCount() const {
        return _propertiesArr[HeavyAtomCount];
    }

    PropertiesFilter::property_t PropertiesFilter::Properties::pubchemCactvsTautoCount() const {
        return _propertiesArr[CactvsTautoCount];
    }

    PropertiesFilter::property_t PropertiesFilter::Properties::pubchemIsotopicAtomCount() const {
        return _propertiesArr[IsotopicAtomCount];
    }

    PropertiesFilter::property_t PropertiesFilter::Properties::pubchemCactvsHbondDonor() const {
        return _propertiesArr[CactvsHbondDonor];
    }

    PropertiesFilter::property_t PropertiesFilter::Properties::pubchemCactvsRotatableBond() const {
        return _propertiesArr[CactvsRotatableBond];
    }

    PropertiesFilter::property_t PropertiesFilter::Properties::pubchemMonoisotopicWeight() const {
        return _propertiesArr[MonoisotopicWeight];
    }

    PropertiesFilter::property_t PropertiesFilter::Properties::pubchemCactvsHbondAcceptor() const {
        return _propertiesArr[CactvsHbondAcceptor];
    }

    PropertiesFilter::property_t PropertiesFilter::Properties::pubchemAtomDefStereoCount() const {
        return _propertiesArr[AtomDefStereoCount];
    }

    PropertiesFilter::property_t PropertiesFilter::Properties::pubchemCompoundCid() const {
        return _propertiesArr[CompoundCid];
    }

    PropertiesFilter::property_t PropertiesFilter::Properties::pubchemMolecularWeight() const {
        return _propertiesArr[MolecularWeight];
    }

    PropertiesFilter::property_t PropertiesFilter::Properties::pubchemBondDefStereoCount() const {
        return _propertiesArr[BondDefStereoCount];
    }

    PropertiesFilter::property_t PropertiesFilter::Properties::pubchemTotalCharge() const {
        return _propertiesArr[TotalCharge];
    }

    PropertiesFilter::property_t PropertiesFilter::Properties::pubchemExactMass() const {
        return _propertiesArr[ExactMass];
    }

    PropertiesFilter::property_t PropertiesFilter::Properties::pubchemCactvsComplexity() const {
        return _propertiesArr[CactvsComplexity];
    }

    PropertiesFilter::property_t PropertiesFilter::Properties::pubchemBondUdefStereoCount() const {
        return _propertiesArr[BondUdefStereoCount];
    }

    PropertiesFilter::property_t PropertiesFilter::Properties::pubchemCactvsTpsa() const {
        return _propertiesArr[CactvsTpsa];
    }

    PropertiesFilter::property_t PropertiesFilter::Properties::pubchemCompoundCanonicalized() const {
        return _propertiesArr[CompoundCanonicalized];
    }

} // qtr