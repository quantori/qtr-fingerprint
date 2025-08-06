#pragma once

#include <memory>
#include <unordered_map>
#include <mutex>
#include <atomic>
#include <chrono>
#include <thread>

#include "search/engines/SearchEngineInterface.h"
#include "search/engines/RDKitSearchEngine.h"
#include "search/engines/IndigoSearchEngine.h"
#include "search/engines/BallTreeSearchEngine.h"
#include "search/engines/SearchEngineType.h"
#include "search/utils/SearchQuery.h"
#include "search/utils/SearchResult.h"
#include "frameworks/RDKitFramework.h"
#include "frameworks/IndigoFramework.h"
#include "dataset/SmilesStorage.h"
#include "utils/Config.h"
#include "CollectionManager.h"

struct TimedSearchResult {
    std::unique_ptr<SearchResult<size_t>> data;
    std::chrono::steady_clock::time_point createdAt;
    bool completed;
    std::string error;
    std::string collectionId;
    
    explicit TimedSearchResult(const std::string& collId = "", std::unique_ptr<SearchResult<size_t>> result = nullptr)
        : data(std::move(result)), createdAt(std::chrono::steady_clock::now()), completed(false), collectionId(collId) {}
};

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
class ServiceManager {
public:
    explicit ServiceManager(Config config);
    
    // Collection management
    std::string createCollection();
    bool addMolecule(const std::string& collectionId, const std::string& smiles);
    bool addMolecules(const std::string& collectionId, const std::vector<std::string>& smilesList);
    bool buildIndex(const std::string& collectionId);
    
    // Collection info
    bool collectionExists(const std::string& collectionId) const;
    CollectionStatus getCollectionStatus(const std::string& collectionId) const;
    std::vector<std::string> getCollectionIds() const;
    size_t getCollectionSize(const std::string& collectionId) const;
    std::string getCollectionError(const std::string& collectionId) const;
    
    // Submit search query and return query ID (collection-specific)
    size_t submitSearch(const std::string& collectionId, const std::string& smiles, size_t maxResults, double timeLimit);
    
    // Get search results by query ID
    std::shared_ptr<TimedSearchResult> getSearchResult(size_t queryId);
    
    // Quick search - immediate response with limited results (collection-specific)
    std::unique_ptr<SearchResult<size_t>> quickSearch(const std::string& collectionId, const std::string& smiles, size_t limit);
    
    // Check if collection is ready for search
    bool isCollectionReady(const std::string& collectionId) const;
    
    // Get SMILES string by molecule ID (collection-specific)
    std::string getSmilesById(const std::string& collectionId, size_t moleculeId) const;
    
    // Cleanup expired queries and old collections
    void cleanupExpiredQueries();
    void cleanupOldCollections();

private:
    Config _config;
    
    // Collection management
    std::unique_ptr<CollectionManager<SearchEngineT>> _collectionManager;
    
    // Query management
    std::mutex _queryMutex;
    std::unordered_map<size_t, std::shared_ptr<TimedSearchResult>> _activeQueries;
    std::atomic<size_t> _queryCounter{0};
    
    // Helper method for executing search
    std::unique_ptr<SearchResult<size_t>> executeSearch(
        const std::string& collectionId, const std::string& smiles, size_t maxResults, bool& stopFlag);
};

// Template implementation
template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
ServiceManager<SearchEngineT>::ServiceManager(Config config)
    : _config(std::move(config)), _collectionManager(std::make_unique<CollectionManager<SearchEngineT>>(config)) {
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
std::string ServiceManager<SearchEngineT>::createCollection() {
    return _collectionManager->createCollection();
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
bool ServiceManager<SearchEngineT>::addMolecule(const std::string& collectionId, const std::string& smiles) {
    return _collectionManager->addMolecule(collectionId, smiles);
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
bool ServiceManager<SearchEngineT>::addMolecules(const std::string& collectionId, const std::vector<std::string>& smilesList) {
    return _collectionManager->addMolecules(collectionId, smilesList);
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
bool ServiceManager<SearchEngineT>::buildIndex(const std::string& collectionId) {
    return _collectionManager->buildIndex(collectionId);
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
bool ServiceManager<SearchEngineT>::collectionExists(const std::string& collectionId) const {
    return _collectionManager->collectionExists(collectionId);
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
CollectionStatus ServiceManager<SearchEngineT>::getCollectionStatus(const std::string& collectionId) const {
    return _collectionManager->getCollectionStatus(collectionId);
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
std::vector<std::string> ServiceManager<SearchEngineT>::getCollectionIds() const {
    return _collectionManager->getCollectionIds();
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
size_t ServiceManager<SearchEngineT>::getCollectionSize(const std::string& collectionId) const {
    return _collectionManager->getCollectionSize(collectionId);
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
std::string ServiceManager<SearchEngineT>::getCollectionError(const std::string& collectionId) const {
    return _collectionManager->getCollectionError(collectionId);
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
size_t ServiceManager<SearchEngineT>::submitSearch(const std::string& collectionId, const std::string& smiles, size_t maxResults, double timeLimit) {
    if (!_collectionManager->isCollectionReady(collectionId)) {
        throw std::runtime_error("Collection not ready: " + collectionId);
    }
    
    auto queryId = _queryCounter++;
    auto timedResult = std::make_shared<TimedSearchResult>(collectionId);
    
    {
        std::lock_guard<std::mutex> lock(_queryMutex);
        _activeQueries[queryId] = timedResult;
    }
    
    // Execute search asynchronously
    std::thread([this, queryId, collectionId, smiles, maxResults, timeLimit, timedResult]() {
        try {
            bool stopFlag = false;
            auto result = executeSearch(collectionId, smiles, maxResults, stopFlag);
            
            {
                std::lock_guard<std::mutex> lock(_queryMutex);
                timedResult->data = std::move(result);
                timedResult->completed = true;
            }
            
            LOG(INFO) << "Search completed for query ID: " << queryId << " in collection: " << collectionId;
            
        } catch (const std::exception& e) {
            std::lock_guard<std::mutex> lock(_queryMutex);
            timedResult->error = e.what();
            timedResult->completed = true;
            LOG(ERROR) << "Search failed for query ID " << queryId << " in collection " << collectionId << ": " << e.what();
        }
    }).detach();
    
    return queryId;
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
std::shared_ptr<TimedSearchResult> ServiceManager<SearchEngineT>::getSearchResult(size_t queryId) {
    std::lock_guard<std::mutex> lock(_queryMutex);
    auto it = _activeQueries.find(queryId);
    if (it != _activeQueries.end()) {
        return it->second;
    }
    return nullptr;
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
std::unique_ptr<SearchResult<size_t>> ServiceManager<SearchEngineT>::quickSearch(const std::string& collectionId, const std::string& smiles, size_t limit) {
    if (!_collectionManager->isCollectionReady(collectionId)) {
        throw std::runtime_error("Collection not ready: " + collectionId);
    }
    
    bool stopFlag = false;
    return executeSearch(collectionId, smiles, limit, stopFlag);
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
bool ServiceManager<SearchEngineT>::isCollectionReady(const std::string& collectionId) const {
    return _collectionManager->isCollectionReady(collectionId);
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
std::string ServiceManager<SearchEngineT>::getSmilesById(const std::string& collectionId, size_t moleculeId) const {
    return _collectionManager->getSmilesById(collectionId, moleculeId);
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
void ServiceManager<SearchEngineT>::cleanupExpiredQueries() {
    const auto ttl = std::chrono::minutes(10);
    auto now = std::chrono::steady_clock::now();
    
    std::lock_guard<std::mutex> lock(_queryMutex);
    for (auto it = _activeQueries.begin(); it != _activeQueries.end(); ) {
        if (now - it->second->createdAt > ttl) {
            LOG(INFO) << "Cleaning up expired query ID: " << it->first;
            it = _activeQueries.erase(it);
        } else {
            ++it;
        }
    }
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
void ServiceManager<SearchEngineT>::cleanupOldCollections() {
    _collectionManager->cleanupOldCollections();
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
std::unique_ptr<SearchResult<size_t>> ServiceManager<SearchEngineT>::executeSearch(
    const std::string& collectionId, const std::string& smiles, size_t maxResults, bool& stopFlag) {
    
    // Get the search engine for the collection
    auto* searchEngine = _collectionManager->getSearchEngine(collectionId);
    if (!searchEngine) {
        throw std::runtime_error("No search engine available for collection: " + collectionId);
    }
    
    SearchQuery query(smiles, maxResults, stopFlag);
    auto result = searchEngine->search(query);
    return std::move(result);
} 