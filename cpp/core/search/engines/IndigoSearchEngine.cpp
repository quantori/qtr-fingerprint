#include "IndigoSearchEngine.h"
#include <tbb/parallel_for.h>

#include <random>
#include <execution>
#include <ranges>
#include <cassert>
#include <utility>

#include <glog/logging.h>

#include "indigo.h"
#include "indigo_internal.h"
#include "IndigoException.h"

namespace {
    std::mt19937 random_generator(0);

    std::filesystem::path generateDBPath() {
        std::filesystem::path res;
        do {
            uint32_t dirId = random_generator();
            std::string dirName = "bingoDb_" + std::to_string(dirId);
            res = std::filesystem::temp_directory_path() / dirName;
        } while (std::filesystem::exists(res));
//        LOG(INFO) << "Generated name: " << res;
        return res;
    }
}

IndigoSearchEngine::IndigoSearchEngine(IndigoSearchEngine::FrameworkT framework, SmilesStorage &&dataset,
                                       const Config &config) : IndigoSearchEngine(std::move(framework)) {
    tbb::parallel_for(
            tbb::blocked_range<size_t>(0, dataset.size()),
            [&](const tbb::blocked_range<size_t> &r) {
                for (size_t idx = r.begin(); idx != r.end(); ++idx) {
                    const auto &smiles = dataset.smiles(idx);
                    try {
                        auto mol = _framework.moleculeFromSmiles(smiles);
                        _db.insertRecord(*mol);
                    }
                    catch (const indigo_cpp::IndigoException &e) {
                        LOG(ERROR) << "Error processing smiles " << smiles << ": " << e.what();
                    }
                }
            }
    );
}

IndigoSearchEngine::IndigoSearchEngine(IndigoSearchEngine::FrameworkT framework) : _framework(framework),
                                                                                   _dbFilePath(generateDBPath()), _db(
                indigo_cpp::BingoMolecule::createDatabaseFile(framework.getSession(), _dbFilePath, "")) {
}

std::unique_ptr<SearchResult<IndigoSearchEngine::ResultT>> IndigoSearchEngine::search(const SearchQuery &query) const {
    auto result = std::make_unique<SearchResult<ResultT>>();
    auto queryMol = _framework.queryMoleculeFromSmiles(query.smiles());
    auto subMatcher = _db.searchSub(*queryMol, "");
    for (auto &searchResult: subMatcher) {
        if (query.checkStopFlag()) {
            break;
        }
        // TODO: fix abstract IDs once Bingo NoSQL is fixed: https://github.com/epam/Indigo/issues/2864
        //        result->addResult(mol.getTarget());
        result->addResult((size_t)searchResult.getId());
        if (result->size() >= query.maxResults()) {
            break;
        }
    }
    return result;
}

IndigoSearchEngine::~IndigoSearchEngine() {
    std::filesystem::remove_all(_dbFilePath);
}

StatTable IndigoSearchEngine::getStat() const {
    // Statistics for IndigoSearchEngine is not collected
    return StatTable();
}

std::string IndigoSearchEngine::resultToSmiles(const ResultT &result) const {
    throw std::runtime_error("IndigoSearchEngine::resultToSmiles is not implemented");
    return std::string();
}
