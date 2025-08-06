#pragma once

#include <memory>
#include <unordered_map>
#include <mutex>
#include <shared_mutex>
#include <atomic>
#include <string>
#include <vector>
#include <chrono>
#include <thread>
#include <glog/logging.h>

#include "search/engines/SearchEngineInterface.h"
#include "dataset/SmilesStorage.h"
#include "utils/Config.h"

enum class CollectionStatus {
    Building,    // Collection exists but molecules can still be added
    Ready,       // Index built, ready for search
    Building_Index, // Currently building index (long running operation)
    Error        // Error state
};

struct Collection {
    std::string id;
    std::vector<std::string> molecules;
    CollectionStatus status;
    std::string errorMessage;
    std::chrono::steady_clock::time_point createdAt;
    std::chrono::steady_clock::time_point lastModified;
    
    Collection(const std::string& collectionId) 
        : id(collectionId), status(CollectionStatus::Building), 
          createdAt(std::chrono::steady_clock::now()),
          lastModified(std::chrono::steady_clock::now()) {}
};

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
class CollectionManager {
public:
    explicit CollectionManager(Config config);
    
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
    
    // Search operations
    SearchEngineT* getSearchEngine(const std::string& collectionId);
    bool isCollectionReady(const std::string& collectionId) const;
    
    // Cleanup
    void cleanupOldCollections();
    
    // Get SMILES by ID within collection
    std::string getSmilesById(const std::string& collectionId, size_t moleculeId) const;

private:
    Config _config;
    mutable std::shared_mutex _collectionsMutex;
    std::unordered_map<std::string, std::unique_ptr<Collection>> _collections;
    std::unordered_map<std::string, std::unique_ptr<SearchEngineT>> _searchEngines;
    std::atomic<size_t> _collectionCounter{0};
    
    // Helper methods
    std::string generateCollectionId();
    void buildIndexAsync(const std::string& collectionId);
};

// Template implementation
template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
CollectionManager<SearchEngineT>::CollectionManager(Config config)
    : _config(std::move(config)) {
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
std::string CollectionManager<SearchEngineT>::createCollection() {
    auto collectionId = generateCollectionId();
    auto collection = std::make_unique<Collection>(collectionId);
    
    {
        std::unique_lock<std::shared_mutex> lock(_collectionsMutex);
        _collections[collectionId] = std::move(collection);
    }
    
    LOG(INFO) << "Created collection: " << collectionId;
    return collectionId;
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
bool CollectionManager<SearchEngineT>::addMolecule(const std::string& collectionId, const std::string& smiles) {
    std::unique_lock<std::shared_mutex> lock(_collectionsMutex);
    
    auto it = _collections.find(collectionId);
    if (it == _collections.end()) {
        return false;
    }
    
    auto& collection = it->second;
    if (collection->status != CollectionStatus::Building) {
        return false; // Can only add molecules to building collections
    }
    
    collection->molecules.push_back(smiles);
    collection->lastModified = std::chrono::steady_clock::now();
    
    LOG(INFO) << "Added molecule to collection " << collectionId << ". Total molecules: " << collection->molecules.size();
    return true;
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
bool CollectionManager<SearchEngineT>::addMolecules(const std::string& collectionId, const std::vector<std::string>& smilesList) {
    std::unique_lock<std::shared_mutex> lock(_collectionsMutex);
    
    auto it = _collections.find(collectionId);
    if (it == _collections.end()) {
        return false;
    }
    
    auto& collection = it->second;
    if (collection->status != CollectionStatus::Building) {
        return false; // Can only add molecules to building collections
    }
    
    collection->molecules.insert(collection->molecules.end(), smilesList.begin(), smilesList.end());
    collection->lastModified = std::chrono::steady_clock::now();
    
    LOG(INFO) << "Added " << smilesList.size() << " molecules to collection " << collectionId 
              << ". Total molecules: " << collection->molecules.size();
    return true;
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
bool CollectionManager<SearchEngineT>::buildIndex(const std::string& collectionId) {
    std::unique_lock<std::shared_mutex> lock(_collectionsMutex);
    
    auto it = _collections.find(collectionId);
    if (it == _collections.end()) {
        return false;
    }
    
    auto& collection = it->second;
    if (collection->status != CollectionStatus::Building) {
        return false; // Can only build index for building collections
    }
    
    if (collection->molecules.empty()) {
        collection->status = CollectionStatus::Error;
        collection->errorMessage = "Cannot build index for empty collection";
        return false;
    }
    
    collection->status = CollectionStatus::Building_Index;
    
    // Start async index building
    std::thread([this, collectionId]() {
        buildIndexAsync(collectionId);
    }).detach();
    
    LOG(INFO) << "Started building index for collection " << collectionId 
              << " with " << collection->molecules.size() << " molecules";
    return true;
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
bool CollectionManager<SearchEngineT>::collectionExists(const std::string& collectionId) const {
    std::shared_lock<std::shared_mutex> lock(_collectionsMutex);
    return _collections.find(collectionId) != _collections.end();
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
CollectionStatus CollectionManager<SearchEngineT>::getCollectionStatus(const std::string& collectionId) const {
    std::shared_lock<std::shared_mutex> lock(_collectionsMutex);
    
    auto it = _collections.find(collectionId);
    if (it == _collections.end()) {
        throw std::runtime_error("Collection not found: " + collectionId);
    }
    
    return it->second->status;
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
std::vector<std::string> CollectionManager<SearchEngineT>::getCollectionIds() const {
    std::shared_lock<std::shared_mutex> lock(_collectionsMutex);
    
    std::vector<std::string> ids;
    ids.reserve(_collections.size());
    
    for (const auto& pair : _collections) {
        ids.push_back(pair.first);
    }
    
    return ids;
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
size_t CollectionManager<SearchEngineT>::getCollectionSize(const std::string& collectionId) const {
    std::shared_lock<std::shared_mutex> lock(_collectionsMutex);
    
    auto it = _collections.find(collectionId);
    if (it == _collections.end()) {
        return 0;
    }
    
    return it->second->molecules.size();
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
std::string CollectionManager<SearchEngineT>::getCollectionError(const std::string& collectionId) const {
    std::shared_lock<std::shared_mutex> lock(_collectionsMutex);
    
    auto it = _collections.find(collectionId);
    if (it == _collections.end()) {
        return "Collection not found";
    }
    
    return it->second->errorMessage;
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
SearchEngineT* CollectionManager<SearchEngineT>::getSearchEngine(const std::string& collectionId) {
    std::shared_lock<std::shared_mutex> lock(_collectionsMutex);
    
    auto it = _searchEngines.find(collectionId);
    if (it != _searchEngines.end()) {
        return it->second.get();
    }
    
    return nullptr;
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
bool CollectionManager<SearchEngineT>::isCollectionReady(const std::string& collectionId) const {
    std::shared_lock<std::shared_mutex> lock(_collectionsMutex);
    
    auto it = _collections.find(collectionId);
    if (it == _collections.end()) {
        return false;
    }
    
    return it->second->status == CollectionStatus::Ready;
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
std::string CollectionManager<SearchEngineT>::getSmilesById(const std::string& collectionId, size_t moleculeId) const {
    std::shared_lock<std::shared_mutex> lock(_collectionsMutex);
    
    auto collectionIt = _collections.find(collectionId);
    if (collectionIt == _collections.end()) {
        return "";
    }
    
    const auto& collection = collectionIt->second;
    
    // If collection is still building, use the temporary molecules vector
    if (collection->status == CollectionStatus::Building) {
        if (moleculeId >= collection->molecules.size()) {
            return "";
        }
        return collection->molecules[moleculeId];
    }
    
    // If collection is ready, use the search engine's resultToSmiles method
    auto engineIt = _searchEngines.find(collectionId);
    if (engineIt == _searchEngines.end()) {
        return "";
    }
    
    try {
        // Use the SearchEngine interface to convert molecule ID to SMILES
        return engineIt->second->resultToSmiles(moleculeId);
    } catch (const std::exception& e) {
        LOG(ERROR) << "Failed to get SMILES for molecule ID " << moleculeId 
                   << " in collection " << collectionId << ": " << e.what();
        return "";
    }
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
void CollectionManager<SearchEngineT>::cleanupOldCollections() {
    const auto ttl = std::chrono::hours(24); // Keep collections for 24 hours
    auto now = std::chrono::steady_clock::now();
    
    std::unique_lock<std::shared_mutex> lock(_collectionsMutex);
    
    for (auto it = _collections.begin(); it != _collections.end(); ) {
        if (now - it->second->lastModified > ttl) {
            LOG(INFO) << "Cleaning up old collection: " << it->first;
            _searchEngines.erase(it->first);
            it = _collections.erase(it);
        } else {
            ++it;
        }
    }
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
std::string CollectionManager<SearchEngineT>::generateCollectionId() {
    auto counter = _collectionCounter++;
    auto timestamp = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::system_clock::now().time_since_epoch()).count();
    
    return "collection_" + std::to_string(timestamp) + "_" + std::to_string(counter);
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
void CollectionManager<SearchEngineT>::buildIndexAsync(const std::string& collectionId) {
    try {
        std::vector<std::string> molecules;
        
        // Move molecules under lock to avoid copying
        {
            std::unique_lock<std::shared_mutex> lock(_collectionsMutex);
            auto it = _collections.find(collectionId);
            if (it == _collections.end()) {
                LOG(ERROR) << "Collection not found during index building: " << collectionId;
                return;
            }
            // Move molecules to avoid copying
            molecules = std::move(it->second->molecules);
        }
        
        LOG(INFO) << "Building search engine for collection " << collectionId 
                  << " with " << molecules.size() << " molecules";
        
        // Create SmilesStorage directly from moved vector
        SmilesStorage storage(std::move(molecules));
        
        // Initialize framework
        auto& framework = SearchEngineT::FrameworkT::getInstance();
        framework.init(_config);
        
        // Create search engine with SmilesStorage
        auto searchEngine = std::make_unique<SearchEngineT>(framework, std::move(storage), _config);
        
        // Update collection status under lock
        {
            std::unique_lock<std::shared_mutex> lock(_collectionsMutex);
            auto it = _collections.find(collectionId);
            if (it != _collections.end()) {
                it->second->status = CollectionStatus::Ready;
                it->second->lastModified = std::chrono::steady_clock::now();
                // Clear the temporary molecules vector since they're now in the search engine
                it->second->molecules.clear();
                it->second->molecules.shrink_to_fit();
                
                _searchEngines[collectionId] = std::move(searchEngine);
                LOG(INFO) << "Successfully built index for collection " << collectionId;
            }
        }
        
    } catch (const std::exception& e) {
        LOG(ERROR) << "Failed to build index for collection " << collectionId << ": " << e.what();
        
        std::unique_lock<std::shared_mutex> lock(_collectionsMutex);
        auto it = _collections.find(collectionId);
        if (it != _collections.end()) {
            it->second->status = CollectionStatus::Error;
            it->second->errorMessage = e.what();
        }
    }
}