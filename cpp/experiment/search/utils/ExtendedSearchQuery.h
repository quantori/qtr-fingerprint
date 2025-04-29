#pragma once

#include "search/utils/SearchQuery.h"
#include "frameworks/FrameworkInterface.h"

template<typename FrameworkT> requires FrameworkInterface<FrameworkT>
class ExtendedSearchQuery : public SearchQuery {
public:
    explicit ExtendedSearchQuery(SearchQuery query) : SearchQuery(query) {
        _queryMolecule = FrameworkT::queryMoleculeFromSmiles(_smiles);
        auto molecule = FrameworkT::moleculeFromSmiles(_smiles);
        _queryFingerprint = FrameworkT::fingerprintFromMolecule(*molecule);
    }

    const typename FrameworkT::QueryMoleculeT &molecule() const {
        return *_queryMolecule;
    }

    const typename FrameworkT::FingerprintT &fingerprint() const {
        return *_queryFingerprint;
    }

private:
    std::unique_ptr<typename FrameworkT::QueryMoleculeT> _queryMolecule;
    std::unique_ptr<typename FrameworkT::FingerprintT> _queryFingerprint;
};
