#include "RDKitSearchEngine.h"

#include <execution>
#include <ranges>
#include <algorithm>

#include <glog/logging.h>

#include "search/utils/utils.h"

namespace {
    void addMol(std::unique_ptr<RDKitFramework::MoleculeT> &&mol,
                const boost::shared_ptr<RDKit::CachedMolHolder> &molHandler,
                const boost::shared_ptr<RDKit::PatternHolder> &fpHandler,
                std::mutex &mutex
    ) {
        auto fp = RDKitFramework::fingerprintFromMolecule(*mol);
        auto compressedMol = RDKitFramework::compressMolecule(*mol);
        {
            std::lock_guard<std::mutex> lockGuard(mutex);
            fpHandler->addFingerprint(*fp);
            molHandler->addBinary(*compressedMol);
        }
    }
}

std::unique_ptr<RDKitSearchEngine::ResultT> RDKitSearchEngine::search(const SearchQuery &query) const {
    auto queryMol = FrameworkT::queryMoleculeFromSmiles(query.smiles());
    auto fp = FrameworkT::fingerprintFromMolecule(*queryMol);
    auto result = std::make_unique<ResultT>();
    int maxResults = query.maxResults() == std::numeric_limits<size_t>::max() ? -1 : (int) query.maxResults();
    auto params = RDKitFramework::getSubstructMatchParameters();
    const unsigned int SearchBlockSize = 100000;

    for (unsigned int block = 0; block < _substructLibrary->size(); block += SearchBlockSize) {
        if (checkShouldStopSearch(query, *result)) {
            break;
        }
        auto matches = _substructLibrary->getMatches(*queryMol, block, block + SearchBlockSize, params, 1, maxResults);
        maxResults -= (int) matches.size();
        for (size_t matchIdx: matches) {
            result->addResult(*_substructLibrary->getMol(matchIdx));
        }
    }
    return result;
}

RDKitSearchEngine::RDKitSearchEngine(SmilesStorage &&dataset) {
    auto molHandler = boost::make_shared<RDKit::CachedMolHolder>();
    auto fpHandler = boost::make_shared<RDKit::PatternHolder>();
    std::mutex mutex;
    auto range = std::views::iota(size_t(0), dataset.size());
    std::for_each(std::execution::par, range.begin(), range.end(), [&](size_t idx) {
        auto &s = dataset.smiles(idx);
        auto mol = RDKitFramework::moleculeFromSmiles(s);
        if (mol == nullptr) {
            LOG(WARNING) << "Can't parse smiles: " << s;
            return;
        }
        addMol(std::move(mol), molHandler, fpHandler, mutex);
    });
    _substructLibrary = std::make_unique<RDKit::SubstructLibrary>(molHandler, fpHandler);
}

