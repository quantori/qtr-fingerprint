#include "RDKitSearchEngine.h"

#include <execution>
#include <ranges>
#include <algorithm>

#include <glog/logging.h>

namespace {
    void addMol(std::unique_ptr<RDKitFramework::MoleculeT> &&mol,
                const boost::shared_ptr<RDKit::CachedMolHolder> &molHandler,
                const boost::shared_ptr<RDKit::PatternHolder> &fpHandler,
                std::mutex &mutex
    ) {
        auto fp = std::unique_ptr<ExplicitBitVect>(
                RDKit::PatternFingerprintMol(*mol, 2048) // Cannot create custom fingerprints
        );
        auto compressedMol = RDKitFramework::compressMolecule(*mol);
        {
            std::lock_guard<std::mutex> lockGuard(mutex);
            fpHandler->addFingerprint(*fp);
            molHandler->addBinary(*compressedMol);
        }
    }
}

std::unique_ptr<SearchResult<RDKitSearchEngine::ResultT>> RDKitSearchEngine::search(const SearchQuery &query) const {
    auto queryMol = FrameworkT::queryMoleculeFromSmiles(query.smiles());
    auto result = std::make_unique<SearchResult<ResultT>>();
    int maxResults = query.maxResults() == std::numeric_limits<size_t>::max() ? -1 : (int) query.maxResults();
    auto params = RDKitFramework::getSubstructMatchParameters();
    const unsigned int SearchBlockSize = 1000000;
    const bool &stopFlag = query.stopFlag();
    for (unsigned int block = 0;
         block < _substructLibrary->size() && !stopFlag && maxResults > 0; block += SearchBlockSize) {
        auto matches = _substructLibrary->getMatches(*queryMol, block, block + SearchBlockSize, params, 1, maxResults);
        maxResults -= (int) matches.size();
        for (size_t matchIdx: matches) {
            result->addResult(matchIdx);
        }
    }
    return result;
}

RDKitSearchEngine::RDKitSearchEngine(RDKitSearchEngine::FrameworkT framework, SmilesStorage &&dataset,
                                     const Config &config) : _framework(framework) {
    auto molHandler = boost::make_shared<RDKit::CachedMolHolder>();
    auto fpHandler = boost::make_shared<RDKit::PatternHolder>();
    std::mutex mutex;
    auto range = std::views::iota(size_t(0), dataset.size());
    std::for_each(std::execution::par, range.begin(), range.end(), [&](size_t idx) {
        auto &s = dataset.smiles(idx);
        auto mol = _framework.moleculeFromSmiles(s);
        if (mol == nullptr) {
            LOG(WARNING) << "Can't parse smiles: " << s;
            return;
        }
        addMol(std::move(mol), molHandler, fpHandler, mutex);
    });
    _substructLibrary = std::make_unique<RDKit::SubstructLibrary>(molHandler, fpHandler);
}

StatTable RDKitSearchEngine::getStat() {
    // Statistics for RDKitSearchEngine is not collected
    return {};
}

std::string RDKitSearchEngine::resultToSmiles(const ResultT &result) const {
    auto mol = _substructLibrary->getMol(result);
    auto smiles = _framework.moleculeToSmiles(*mol);
    return smiles;
}

