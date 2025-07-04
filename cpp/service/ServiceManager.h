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

struct TimedSearchResult {
    std::unique_ptr<SearchResult<size_t>> data;
    std::chrono::steady_clock::time_point createdAt;
    bool completed;
    std::string error;
    
    explicit TimedSearchResult(std::unique_ptr<SearchResult<size_t>> result = nullptr)
        : data(std::move(result)), createdAt(std::chrono::steady_clock::now()), completed(false) {}
};

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
class ServiceManager {
public:
    explicit ServiceManager(Config config);
    
    // Initialize search engine with dataset
    void initializeEngine(SmilesStorage&& dataset);
    
    // Submit search query and return query ID
    size_t submitSearch(const std::string& smiles, size_t maxResults, double timeLimit);
    
    // Get search results by query ID
    std::shared_ptr<TimedSearchResult> getSearchResult(size_t queryId);
    
    // Quick search - immediate response with limited results
    std::unique_ptr<SearchResult<size_t>> quickSearch(const std::string& smiles, size_t limit);
    
    // Check if engine is ready
    bool isEngineReady() const;
    
    // Convert index to SMILES string
    std::string indexToString(size_t index) const;
    
    // Get SMILES string by molecule ID
    std::string getSmilesById(size_t moleculeId) const;
    
    // Cleanup expired queries
    void cleanupExpiredQueries();

private:
    Config _config;
    std::atomic<bool> _engineReady{false};
    
    // Query management
    std::mutex _queryMutex;
    std::unordered_map<size_t, std::shared_ptr<TimedSearchResult>> _activeQueries;
    std::atomic<size_t> _queryCounter{0};
    
    // Single templated search engine
    std::unique_ptr<SearchEngineT> _searchEngine;
    
    // Helper method for executing search
    std::unique_ptr<SearchResult<size_t>> executeSearch(
        const std::string& smiles, size_t maxResults, bool& stopFlag);
};

// Template implementation
template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
ServiceManager<SearchEngineT>::ServiceManager(Config config)
    : _config(std::move(config)) {
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
void ServiceManager<SearchEngineT>::initializeEngine(SmilesStorage&& dataset) {
    LOG(INFO) << "Initializing search engine";
    
    _engineReady = false;
    _searchEngine.reset();
    
    try {
        auto& framework = SearchEngineT::FrameworkT::getInstance();
        framework.init(_config);
        
        _searchEngine = std::make_unique<SearchEngineT>(framework, std::move(dataset), _config);
        
        _engineReady = true;
        LOG(INFO) << "Search engine initialized successfully";
        
    } catch (const std::exception& e) {
        LOG(ERROR) << "Failed to initialize search engine: " << e.what();
        throw;
    }
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
size_t ServiceManager<SearchEngineT>::submitSearch(const std::string& smiles, size_t maxResults, double timeLimit) {
    if (!_engineReady) {
        throw std::runtime_error("Search engine not ready");
    }
    
    auto queryId = _queryCounter++;
    auto timedResult = std::make_shared<TimedSearchResult>();
    
    {
        std::lock_guard<std::mutex> lock(_queryMutex);
        _activeQueries[queryId] = timedResult;
    }
    
    // Execute search asynchronously
    std::thread([this, queryId, smiles, maxResults, timeLimit, timedResult]() {
        try {
            bool stopFlag = false;
            auto result = executeSearch(smiles, maxResults, stopFlag);
            
            {
                std::lock_guard<std::mutex> lock(_queryMutex);
                timedResult->data = std::move(result);
                timedResult->completed = true;
            }
            
            LOG(INFO) << "Search completed for query ID: " << queryId;
            
        } catch (const std::exception& e) {
            std::lock_guard<std::mutex> lock(_queryMutex);
            timedResult->error = e.what();
            timedResult->completed = true;
            LOG(ERROR) << "Search failed for query ID " << queryId << ": " << e.what();
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
std::unique_ptr<SearchResult<size_t>> ServiceManager<SearchEngineT>::quickSearch(const std::string& smiles, size_t limit) {
    if (!_engineReady) {
        throw std::runtime_error("Search engine not ready");
    }
    
    bool stopFlag = false;
    return executeSearch(smiles, limit, stopFlag);
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
bool ServiceManager<SearchEngineT>::isEngineReady() const {
    return _engineReady;
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
std::string ServiceManager<SearchEngineT>::indexToString(size_t index) const {
    // TODO: Implement actual index to SMILES conversion
    return std::string();
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
std::string ServiceManager<SearchEngineT>::getSmilesById(size_t moleculeId) const {
    if (!_engineReady || !_searchEngine) {
        throw std::runtime_error("Search engine not ready");
    }
    
    try {
        // Use the search engine's resultToSmiles method to convert molecule ID to SMILES
        return _searchEngine->resultToSmiles(moleculeId);
    } catch (const std::exception& e) {
        LOG(ERROR) << "Failed to get SMILES for molecule ID " << moleculeId << ": " << e.what();
        return std::string(); // Return empty string if not found or error
    }
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
std::unique_ptr<SearchResult<size_t>> ServiceManager<SearchEngineT>::executeSearch(
    const std::string& smiles, size_t maxResults, bool& stopFlag) {
    
    SearchQuery query(smiles, maxResults, stopFlag);
    auto result = _searchEngine->search(query);
    return std::move(result);
} 