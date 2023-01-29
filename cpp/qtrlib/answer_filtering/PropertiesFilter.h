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

        struct Properties {
            property_t pubchemComponentCount;
            property_t pubchemXlogp3;
            property_t pubchemAtomUdefStereoCount;
            property_t pubchemHeavyAtomCount;
            property_t pubchemCactvsTautoCount;
            property_t pubchemIsotopicAtomCount;
            property_t pubchemCactvsHbondDonor;
            property_t pubchemCactvsRotatableBond;
            property_t pubchemMonoisotopicWeight;
            property_t pubchemCactvsHbondAcceptor;
            property_t pubchemAtomDefStereoCount;
            property_t pubchemCompoundCid;
            property_t pubchemMolecularWeight;
            property_t pubchemBondDefStereoCount;
            property_t pubchemTotalCharge;
            property_t pubchemExactMass;
            property_t pubchemCactvsComplexity;
            property_t pubchemBondUdefStereoCount;
            property_t pubchemCactvsTpsa;
            property_t pubchemCompoundCanonicalized;

            Properties() = default;

            Properties(const Properties &) = default;

            property_t operator[](size_t i) const;

            property_t &operator[](size_t i);

            static constexpr size_t size();
        };

        struct Bounds {
            Bounds();

            Bounds(const Bounds &) = default;

            [[nodiscard]] bool Check(const Properties &properties) const;

            Properties minBounds{};
            Properties maxBounds{};
        };

        bool operator()(CIDType id) override;

        std::unique_ptr<AnswerFilter> copy() override;

        explicit PropertiesFilter(std::shared_ptr<const std::vector<Properties>> propertiesTable);

        PropertiesFilter(std::shared_ptr<const std::vector<Properties>> propertiesTable, Bounds bounds);

        void setBounds(Bounds bounds);

        void setBounds(std::shared_ptr<const Bounds> bounds);

    private:
        std::shared_ptr<const Bounds> _bounds;
        std::shared_ptr<const std::vector<Properties>> _propertiesTable;
    };


} // qtr
