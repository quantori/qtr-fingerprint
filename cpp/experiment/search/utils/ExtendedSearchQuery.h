#pragma once

#include "search/utils/SearchQuery.h"
#include "frameworks/FrameworkInterface.h"

template<typename FrameworkT> requires FrameworkInterface<FrameworkT>
class ExtendedSearchQuery : public SearchQuery {
public:
    explicit ExtendedSearchQuery(SearchQuery query) : SearchQuery(query) {
        auto &framework = FrameworkT::getInstance();
        _queryMolecule = framework.queryMoleculeFromSmiles(_smiles);
        auto molecule = framework.moleculeFromSmiles(_smiles);
        auto fingerprint = framework.fingerprintFromMolecule(*molecule);
        _queryFingerprint = framework.queryFingerprintFromFingerprint(*fingerprint);
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
