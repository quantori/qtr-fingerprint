#include "RDKitSearchEngine.h"

#include "GraphMol/Substruct/SubstructMatch.h"

#include <glog/logging.h>

#include <future>
#include <ranges>

namespace {

    RDKit::SubstructMatchParameters getSubstructMatchParameters() {
        RDKit::SubstructMatchParameters params;
        params.recursionPossible = true;
        params.useChirality = true;
        params.useQueryQueryMatches = false;
        return params;
    }
}

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
std::vector<uint64_t> RDKitSearchEngine::getMatches(const MoleculeType &queryMol, int maxResults, bool &stopFlag) {
    const unsigned int STEP = 100000;
    std::vector<uint64_t> result;
    for (unsigned int block = 0; block < _substructLibrary->size() && maxResults > 0 && !stopFlag; block += STEP) {
        auto matches = _substructLibrary->getMatches(queryMol, block, block + STEP, getSubstructMatchParameters(), 1,
                                                     maxResults);
        maxResults -= (int) matches.size();
        result.insert(result.end(), matches.begin(), matches.end());
    }
    return result;
}

std::vector<uint64_t> RDKitSearchEngine::getMatches(const RDKitSearchEngine::QueryMoleculeType &queryMol,
                                                    const RDKitSearchEngine::FingerprintType &fingerprint,
                                                    int maxResults, bool &stopFlag) {
    std::vector<uint64_t> result;
    auto &molecules = _substructLibrary->getMolecules();
    auto &fingerprints = _substructLibrary->getFingerprints();
    auto params = getSubstructMatchParameters();
    assert(molecules.size() == fingerprints.size());
    for (uint64_t i = 0; i < molecules.size() && !stopFlag && result.size() < maxResults; i++) {
        if (!fingerprints.passesFilter(i, fingerprint.bitVector())) {
            continue;
        }
        auto mol = molecules.getMol(i);
        if (!SubstructMatch(*mol, queryMol, params).empty()) {
            result.push_back(i);
        }
    }
    return result;
}

std::unique_ptr<RDKitSearchEngine::MoleculeType> RDKitSearchEngine::smilesToMolecule(const std::string &smiles) {
    return std::unique_ptr<MoleculeType>(RDKit::SmilesToMol(smiles));
}

RDKitSearchEngine::RDKitSearchEngine(const std::vector<std::string> &smiles) {
    auto molHandler = boost::make_shared<RDKit::CachedMolHolder>();
    auto fpHandler = boost::make_shared<RDKit::PatternHolder>();
    _substructLibrary = std::make_shared<RDKit::SubstructLibrary>(molHandler, fpHandler);
    for (auto &s: smiles) {
        std::unique_ptr<MoleculeType> mol;
        try {
            mol = smilesToMolecule(s);
        } catch (const std::exception &e) {
            LOG(WARNING) << "Skip smiles: " << s << " error: " << e.what();
            continue;
        }
        if (mol == nullptr) {
            LOG(WARNING) << "Can't parse smiles: " << s;
            continue;
        }
        auto fp = std::make_unique<RDKitFingerprint>(*mol);
        _substructLibrary->addMol(*mol);
    }
}

RDKitSearchEngine::RDKitSearchEngine(std::vector<std::unique_ptr<StorageMoleculeType>> &&molecules) {
    auto molHandler = boost::make_shared<RDKit::CachedMolHolder>();
    auto fpHandler = boost::make_shared<RDKit::PatternHolder>();
    _substructLibrary = std::make_shared<RDKit::SubstructLibrary>(molHandler, fpHandler);
    for (auto &molPickle: molecules) {
        RDKit::ROMol mol;
        RDKit::MolPickler::molFromPickle(*molPickle, mol);
        _substructLibrary->addMol(mol);
    }

}

RDKitSearchEngine::RDKitSearchEngine(
        std::vector<std::pair<std::unique_ptr<StorageMoleculeType>, std::unique_ptr<FingerprintType>>> &&data) {
    auto molHandler = boost::make_shared<RDKit::CachedMolHolder>();
    auto fpHandler = boost::make_shared<RDKit::PatternHolder>();
    for (auto &[molPickle, fp]: data) {
        molHandler->addBinary(*molPickle);
        fpHandler->addFingerprint(fp->bitVector());
    }
    _substructLibrary = std::make_shared<RDKit::SubstructLibrary>(molHandler, fpHandler);
}

std::unique_ptr<RDKitSearchEngine::MoleculeType> RDKitSearchEngine::smilesToQueryMolecule(const std::string &smiles) {
    return smilesToMolecule(smiles);
}

std::unique_ptr<RDKitSearchEngine::StorageMoleculeType>
RDKitSearchEngine::moleculeToStorageMolecule(const RDKitSearchEngine::MoleculeType &molecule) {
    auto storageMol = std::make_unique<StorageMoleculeType>();
    RDKit::MolPickler::pickleMol(molecule, *storageMol);
    return storageMol;
}

std::unique_ptr<RDKitSearchEngine::MoleculeType>
RDKitSearchEngine::storageMoleculeToMolecule(const RDKitSearchEngine::StorageMoleculeType &storageMolecule) {
    auto mol = std::make_unique<MoleculeType>();
    RDKit::MolPickler::molFromPickle(storageMolecule, *mol);
    return mol;
}

