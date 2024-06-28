#include "RDKitSearchEngine.h"

#include <glog/logging.h>

#include <future>
#include <ranges>

/**
 * @brief Searches for substructure matches in the RDKit substructure library.
 *
 * This function searches for substructure matches in the RDKit substructure library
 * using the provided query SMILES string. The search algorithm is aligned with the
 * official RDKit documentation (http://www.rdkit.org/new_docs/source/rdkit.Chem.rdSubstructLibrary.html),
 * with the modification that it checks the stop flag every 100000 iterations.
 *
 * @param queryMol A string; containing the query SMILES representation. // TODO: change doc
 * @param maxResults The maximum number of results to retrieve.
 * @param stopFlag A boolean flag that can be set to true to stop the search prematurely.
 * @return A vector of uint64_t containing the indices of the matching molecules.
 *
 * @note If the query SMILES cannot be parsed, the function throws a std::runtime_error.
 * @note TODO: Implement a version without time limit checks and perform full tests on a large dataset
 *             to ensure that these checks do not affect performance.
 */
std::vector<uint64_t> RDKitSearchEngine::getMatches(const RDKit::ROMol &queryMol, int maxResults, bool &stopFlag) {
    RDKit::SubstructMatchParameters params;
    params.recursionPossible = true;
    params.useChirality = true;
    params.useQueryQueryMatches = false;
    const unsigned int STEP = 100000;
    std::vector<uint64_t> result;
    for (unsigned int block = 0; block < _substructLibrary->size() && maxResults > 0 && !stopFlag; block += STEP) {
        auto matches = _substructLibrary->getMatches(queryMol, block, block + STEP, params, 1, maxResults);
        maxResults -= (int) matches.size();
        result.insert(result.end(), matches.begin(), matches.end());
    }
    return result;
}

std::unique_ptr<RDKit::ROMol> RDKitSearchEngine::smilesToMolecule(const std::string &smiles) {
    return std::unique_ptr<RDKit::ROMol>(RDKit::SmilesToMol(smiles));
}

RDKitSearchEngine::RDKitSearchEngine(const std::vector<std::string> &smiles) {
    auto molHandler = boost::make_shared<RDKit::CachedMolHolder>();
    auto fpHandler = boost::make_shared<RDKit::PatternHolder>();
    _substructLibrary = std::make_shared<RDKit::SubstructLibrary>(molHandler, fpHandler);
    for (auto &s: smiles) {
        auto mol = smilesToMolecule(s);
        _substructLibrary->addMol(*mol);
    }
}

RDKitSearchEngine::RDKitSearchEngine(std::vector<std::unique_ptr<RDKit::ROMol>> &&molecules) {
    auto molHandler = boost::make_shared<RDKit::CachedMolHolder>();
    auto fpHandler = boost::make_shared<RDKit::PatternHolder>();
    _substructLibrary = std::make_shared<RDKit::SubstructLibrary>(molHandler, fpHandler);
    for (auto &mol: molecules) {
        _substructLibrary->addMol(*mol);
    }

}

RDKitSearchEngine::RDKitSearchEngine(
        std::vector<std::pair<std::unique_ptr<RDKit::ROMol>, std::unique_ptr<FingerprintType>>> &&data) {
    auto molHandler = boost::make_shared<RDKit::CachedMolHolder>();
    auto fpHandler = boost::make_shared<RDKit::PatternHolder>();
    for (auto &[mol, fp]: data) {
        molHandler->addMol(*mol);
        fpHandler->addFingerprint(fp->bitVector());
    }
    _substructLibrary = std::make_shared<RDKit::SubstructLibrary>(molHandler, fpHandler);
}
