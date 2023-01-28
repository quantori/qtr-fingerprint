#pragma once

#include "AnswerFilter.h"

#include <numeric>
#include <vector>
#include <string>

namespace qtr {

    class PropertiesFilter : public AnswerFilter {
    public:
        inline static const std::string propertyNames[] = {
                "PUBCHEM_COMPONENT_COUNT",
                "PUBCHEM_XLOGP3",
                "PUBCHEM_ATOM_UDEF_STEREO_COUNT",
                "PUBCHEM_HEAVY_ATOM_COUNT",
                "PUBCHEM_CACTVS_TAUTO_COUNT",
                "PUBCHEM_ISOTOPIC_ATOM_COUNT",
                "PUBCHEM_CACTVS_HBOND_DONOR",
                "PUBCHEM_CACTVS_ROTATABLE_BOND",
                "PUBCHEM_MONOISOTOPIC_WEIGHT",
                "PUBCHEM_CACTVS_HBOND_ACCEPTOR",
                "PUBCHEM_ATOM_DEF_STEREO_COUNT",
                "PUBCHEM_COMPOUND_CID",
                "PUBCHEM_MOLECULAR_WEIGHT",
                "PUBCHEM_BOND_DEF_STEREO_COUNT",
                "PUBCHEM_TOTAL_CHARGE",
                "PUBCHEM_EXACT_MASS",
                "PUBCHEM_CACTVS_COMPLEXITY",
                "PUBCHEM_BOND_UDEF_STEREO_COUNT",
                "PUBCHEM_CACTVS_TPSA",
                "PUBCHEM_COMPOUND_CANONICALIZED"
        };

        using property_t = float;
        static const size_t propertiesCount = std::size(propertyNames);
        using property_list_t = property_t[propertiesCount];

        class Bounds {
        public:
            Bounds();

            Bounds(const Bounds &) = default;

            [[nodiscard]] bool Check(const property_list_t &properties) const;

        private:
            property_list_t _minBounds{};
            property_list_t _maxBounds{};
        };

        bool operator()(CIDType id) override;

        std::unique_ptr<AnswerFilter> copy() override;

        explicit PropertiesFilter(std::shared_ptr<const std::vector<property_list_t>> propertiesTable);

        void setBounds(Bounds bounds);

        void setBounds(std::shared_ptr<const Bounds> bounds);

    private:
        std::shared_ptr<const Bounds> _bounds;
        std::shared_ptr<const std::vector<property_list_t>> _propertiesTable;
    };


} // qtr
