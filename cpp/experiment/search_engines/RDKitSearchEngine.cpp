#include "RDKitSearchEngine.h"

#include <glog/logging.h>

#include <future>
#include <ranges>

namespace {


    void addFileToHandlers(const std::filesystem::path &smilesFile,
                           const boost::shared_ptr <RDKit::CachedMolHolder> &molHandler,
                           const boost::shared_ptr <RDKit::PatternHolder> &patternHolder,
                           std::mutex &mutex) {
        std::ifstream in(smilesFile);
        std::string smiles;
        while (in.peek() != EOF) {
            std::string id;
            in >> id >> smiles;
            std::string lineEnding;
            std::getline(in, lineEnding);
            if (smiles.empty()) {
                LOG(WARNING) << "Found line with wrong formatting and skipped it: \"" << smiles << lineEnding << "\"";
                continue;
            }
            std::unique_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol(smiles));
            if (mol == nullptr) {
                LOG(ERROR) << "Cannot parse smiles: " << smiles;
                continue;
            }
            std::unique_ptr<ExplicitBitVect> fingerprint(RDKit::PatternFingerprintMol(*mol));
            {
                std::lock_guard<std::mutex> lockGuard(mutex);
                molHandler->addMol(*mol);
                patternHolder->addFingerprint(fingerprint.release());
            }
        }
    }
}

RDKitSearchEngine::RDKitSearchEngine(const std::filesystem::path &datasetDir) {
    auto molHandler = boost::make_shared<RDKit::CachedMolHolder>();
    auto fpHandler = boost::make_shared<RDKit::PatternHolder>();
    std::mutex mutex;
    std::vector<std::future<void>> tasks;
    LOG(INFO) << "Start RDKit search engine construction. Dataset dir: " << datasetDir;
    for (auto &dirEntry: std::filesystem::directory_iterator(datasetDir)) {
        auto &filename = dirEntry.path();
        if (!std::filesystem::is_regular_file(filename)) {
            LOG(WARNING) << "Skipped entity: " << filename;
            continue;
        }
        tasks.push_back(
                async(std::launch::async, addFileToHandlers, filename, std::ref(molHandler), std::ref(fpHandler),
                      std::ref(mutex)));
    }
    for (auto &task: tasks) {
        task.wait();
    }

    _substructLibrary = std::make_shared<RDKit::SubstructLibrary>(molHandler, fpHandler);
}

/**
 * @brief Searches for substructure matches in the RDKit substructure library.
 *
 * This function searches for substructure matches in the RDKit substructure library
 * using the provided query SMILES string. The search algorithm is aligned with the
 * official RDKit documentation (http://www.rdkit.org/new_docs/source/rdkit.Chem.rdSubstructLibrary.html),
 * with the modification that it checks the stop flag every 100000 iterations.
 *
 * @param querySmiles A string; containing the query SMILES representation.
 * @param maxResults The maximum number of results to retrieve.
 * @param stopFlag A boolean flag that can be set to true to stop the search prematurely.
 * @return A vector of uint64_t containing the indices of the matching molecules.
 *
 * @note If the query SMILES cannot be parsed, the function throws a std::runtime_error.
 * @note TODO: Implement a version without time limit checks and perform full tests on a large dataset
 *             to ensure that these checks do not affect performance.
 */
std::vector<uint64_t> RDKitSearchEngine::getMatches(const std::string &querySmiles, int maxResults, bool &stopFlag) {
    std::unique_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol(querySmiles));
    if (mol == nullptr) {
        throw std::runtime_error("Cannot parse query smiles: " + querySmiles);
    }
    RDKit::SubstructMatchParameters params;
    params.recursionPossible = true;
    params.useChirality = true;
    params.useQueryQueryMatches = false;
    const unsigned int STEP = 100000;
    std::vector<uint64_t> result;
    for (unsigned int block = 0; block < _substructLibrary->size() && maxResults > 0 && !stopFlag; block += STEP) {
        auto matches = _substructLibrary->getMatches(*mol, block, block + STEP, params, 1, maxResults);
        maxResults -= (int) matches.size();
        result.insert(result.end(), matches.begin(), matches.end());
    }
    return result;
}
