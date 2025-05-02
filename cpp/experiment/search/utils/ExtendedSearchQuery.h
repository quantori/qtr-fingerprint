#pragma once

#include "search/utils/SearchQuery.h"
#include "frameworks/FrameworkInterface.h"

template<typename FrameworkT> requires FrameworkInterface<FrameworkT>
class ExtendedSearchQuery : public SearchQuery {
public:
    explicit ExtendedSearchQuery(SearchQuery query) : SearchQuery(query) {
        _queryMolecule = FrameworkT::queryMoleculeFromSmiles(_smiles);
        auto molecule = FrameworkT::moleculeFromSmiles(_smiles);
        auto fingerprint = FrameworkT::fingerprintFromMolecule(*molecule);
        _queryFingerprint = FrameworkT::queryFingerprintFromFingerprint(*fingerprint);
    }

    const typename FrameworkT::QueryMoleculeT &molecule() const {
        return *_queryMolecule;
    }

    const typename FrameworkT::QueryFingerprintT &fingerprint() const {
        return *_queryFingerprint;
    }

private:
    std::unique_ptr<typename FrameworkT::QueryMoleculeT> _queryMolecule;
    std::unique_ptr<typename FrameworkT::QueryFingerprintT> _queryFingerprint;
};
