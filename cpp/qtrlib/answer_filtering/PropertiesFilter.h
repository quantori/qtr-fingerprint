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

        enum PropertiesIndex {
            ComponentCount,
            Xlogp3,
            AtomUdefStereoCount,
            HeavyAtomCount,
            CactvsTautoCount,
            IsotopicAtomCount,
            CactvsHbondDonor,
            CactvsRotatableBond,
            MonoisotopicWeight,
            CactvsHbondAcceptor,
            AtomDefStereoCount,
            CompoundCid,
            MolecularWeight,
            BondDefStereoCount,
            TotalCharge,
            ExactMass,
            CactvsComplexity,
            BondUdefStereoCount,
            CactvsTpsa,
            CompoundCanonicalized
        };

        using property_t = float;

        class Properties {
        private:
            property_t _propertiesArr[20] = {};

        public:
            [[nodiscard]] property_t pubchemComponentCount() const;

            [[nodiscard]] property_t pubchemXlogp3() const;

            [[nodiscard]] property_t pubchemAtomUdefStereoCount() const;

            [[nodiscard]] property_t pubchemHeavyAtomCount() const;

            [[nodiscard]] property_t pubchemCactvsTautoCount() const;

            [[nodiscard]] property_t pubchemIsotopicAtomCount() const;

            [[nodiscard]] property_t pubchemCactvsHbondDonor() const;

            [[nodiscard]] property_t pubchemCactvsRotatableBond() const;

            [[nodiscard]] property_t pubchemMonoisotopicWeight() const;

            [[nodiscard]] property_t pubchemCactvsHbondAcceptor() const;

            [[nodiscard]] property_t pubchemAtomDefStereoCount() const;

            [[nodiscard]] property_t pubchemCompoundCid() const;

            [[nodiscard]] property_t pubchemMolecularWeight() const;

            [[nodiscard]] property_t pubchemBondDefStereoCount() const;

            [[nodiscard]] property_t pubchemTotalCharge() const;

            [[nodiscard]] property_t pubchemExactMass() const;

            [[nodiscard]] property_t pubchemCactvsComplexity() const;

            [[nodiscard]] property_t pubchemBondUdefStereoCount() const;

            [[nodiscard]] property_t pubchemCactvsTpsa() const;

            [[nodiscard]] property_t pubchemCompoundCanonicalized() const;

            Properties() = default;

            Properties(const Properties &) = default;

            property_t operator[](size_t i) const;

            property_t &operator[](size_t i);

            [[nodiscard]] constexpr size_t size() const {
                return std::size(_propertiesArr);
            }
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
