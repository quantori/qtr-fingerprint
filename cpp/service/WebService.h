#pragma once

#include <memory>
#include <atomic>
#include <thread>
#include <chrono>
#include <sstream>
#include <fstream>
#include <filesystem>
#include <ctime>

#include <crow.h>
#include <glog/logging.h>

#include "ServiceManager.h"
#include "utils/Config.h"
#include "io/SmilesCSVParser.h"

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
class WebService {
public:
    explicit WebService(std::shared_ptr<ServiceManager<SearchEngineT>> serviceManager, int port = 8080);
    
    // Start the web service
    void run();
    
    // Stop the web service
    void stop();

private:
    std::shared_ptr<ServiceManager<SearchEngineT>> _serviceManager;
    crow::SimpleApp _app;
    int _port;
    std::atomic<bool> _shutdownFlag{false};
    std::unique_ptr<std::thread> _cleanupThread;
    
    // Setup all routes
    void setupRoutes();
    
    // Route handlers
    crow::response handlePostQuery(const crow::request& req);
    crow::response handleGetQuery(const crow::request& req);
    crow::response handleFastQuery(const crow::request& req);
    crow::response handleStatus(const crow::request& req);
    crow::response handleSmilesById(const crow::request& req);
    
    // Collection management handlers
    crow::response handleCreateCollection(const crow::request& req);
    crow::response handleAddMolecule(const crow::request& req, const std::string& collectionId);
    crow::response handleUploadCSV(const crow::request& req, const std::string& collectionId);
    crow::response handleBuildIndex(const crow::request& req, const std::string& collectionId);
    crow::response handleGetCollections(const crow::request& req);
    crow::response handleGetCollectionInfo(const crow::request& req, const std::string& collectionId);
    
    // Helper methods
    crow::json::wvalue resultToJson(const SearchResult<size_t>& result, size_t offset = 0, size_t limit = SIZE_MAX);
    std::pair<size_t, size_t> extractBounds(const crow::json::rvalue& json);
    
    // Cleanup thread for expired queries
    void startCleanupThread();
};

// Template implementation
template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
WebService<SearchEngineT>::WebService(std::shared_ptr<ServiceManager<SearchEngineT>> serviceManager, int port)
    : _serviceManager(std::move(serviceManager)), _port(port) {
    setupRoutes();
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
void WebService<SearchEngineT>::run() {
    LOG(INFO) << "Starting web service on port " << _port;
    
    startCleanupThread();
    
    _app.port(_port).concurrency(4).run();
    
    _shutdownFlag = true;
    if (_cleanupThread && _cleanupThread->joinable()) {
        _cleanupThread->join();
    }
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
void WebService<SearchEngineT>::stop() {
    LOG(INFO) << "Stopping web service";
    _shutdownFlag = true;
    _app.stop();
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
void WebService<SearchEngineT>::setupRoutes() {
    // CORS middleware
    CROW_ROUTE(_app, "/").methods(crow::HTTPMethod::OPTIONS)([](const crow::request&) {
        crow::response res(200);
        res.add_header("Access-Control-Allow-Origin", "*");
        res.add_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS");
        res.add_header("Access-Control-Allow-Headers", "Content-Type, Authorization");
        return res;
    });
    
    // POST /query - Submit search query for async processing
    CROW_ROUTE(_app, "/query").methods(crow::HTTPMethod::POST)([this](const crow::request& req) {
        return handlePostQuery(req);
    });
    
    // GET /query - Get search results by query ID
    CROW_ROUTE(_app, "/query").methods(crow::HTTPMethod::GET)([this](const crow::request& req) {
        return handleGetQuery(req);
    });
    
    // POST /fastquery - Quick search with immediate response
    CROW_ROUTE(_app, "/fastquery").methods(crow::HTTPMethod::POST)([this](const crow::request& req) {
        return handleFastQuery(req);
    });
    
    // GET /status - Service status
    CROW_ROUTE(_app, "/status").methods(crow::HTTPMethod::GET)([this](const crow::request& req) {
        return handleStatus(req);
    });
    
    // GET /smiles - Get SMILES by molecule ID
    CROW_ROUTE(_app, "/smiles").methods(crow::HTTPMethod::GET)([this](const crow::request& req) {
        return handleSmilesById(req);
    });
    
    // Collection management routes
    // POST /collections - Create new collection
    CROW_ROUTE(_app, "/collections").methods(crow::HTTPMethod::POST)([this](const crow::request& req) {
        return handleCreateCollection(req);
    });
    
    // GET /collections - List all collections
    CROW_ROUTE(_app, "/collections").methods(crow::HTTPMethod::GET)([this](const crow::request& req) {
        return handleGetCollections(req);
    });
    
    // POST /collections/<id>/molecules - Add molecule(s) to collection
    CROW_ROUTE(_app, "/collections/<string>/molecules").methods(crow::HTTPMethod::POST)([this](const crow::request& req, const std::string& collectionId) {
        return handleAddMolecule(req, collectionId);
    });
    
    // POST /collections/<id>/upload - Upload CSV file with molecules to collection  
    CROW_ROUTE(_app, "/collections/<string>/upload").methods(crow::HTTPMethod::POST)([this](const crow::request& req, const std::string& collectionId) {
        return handleUploadCSV(req, collectionId);
    });
    
    // POST /collections/<id>/build - Build index for collection
    CROW_ROUTE(_app, "/collections/<string>/build").methods(crow::HTTPMethod::POST)([this](const crow::request& req, const std::string& collectionId) {
        return handleBuildIndex(req, collectionId);
    });
    
    // GET /collections/<id> - Get collection info
    CROW_ROUTE(_app, "/collections/<string>").methods(crow::HTTPMethod::GET)([this](const crow::request& req, const std::string& collectionId) {
        return handleGetCollectionInfo(req, collectionId);
    });
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
crow::response WebService<SearchEngineT>::handlePostQuery(const crow::request& req) {
    crow::response res;
    res.add_header("Access-Control-Allow-Origin", "*");
    
    try {
        auto json = crow::json::load(req.body);
        if (!json || !json.has("smiles") || !json.has("collectionId")) {
            res.code = 400;
            res.body = "Missing or invalid 'smiles' or 'collectionId' field";
            return res;
        }
        
        std::string smiles = json["smiles"].s();
        std::string collectionId = json["collectionId"].s();
        auto maxResults = json.has("maxResults") ? static_cast<size_t>(json["maxResults"].i()) : 1000;
        auto timeLimit = json.has("timeLimit") ? json["timeLimit"].d() : 60.0;
        
        if (!_serviceManager->isCollectionReady(collectionId)) {
            res.code = 503;
            res.body = "Collection not ready: " + collectionId;
            return res;
        }
        
        auto queryId = _serviceManager->submitSearch(collectionId, smiles, maxResults, timeLimit);
        
        res.code = 200;
        res.body = std::to_string(queryId);
        
        LOG(INFO) << "Submitted query with ID: " << queryId << " for SMILES: " << smiles << " in collection: " << collectionId;
        
    } catch (const std::exception& e) {
        LOG(ERROR) << "Error in handlePostQuery: " << e.what();
        res.code = 500;
        res.body = "Internal server error: " + std::string(e.what());
    }
    
    return res;
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
crow::response WebService<SearchEngineT>::handleGetQuery(const crow::request& req) {
    crow::response res;
    res.add_header("Access-Control-Allow-Origin", "*");
    
    try {
        const char* idStr = req.url_params.get("searchId");
        const char* offsetStr = req.url_params.get("offset");
        const char* limitStr = req.url_params.get("limit");
        
        if (!idStr) {
            res.code = 400;
            res.body = "Missing 'searchId' parameter";
            return res;
        }
        
        size_t searchId = std::stoull(idStr);
        size_t offset = offsetStr ? std::stoull(offsetStr) : 0;
        size_t limit = limitStr ? std::stoull(limitStr) : 100;
        
        auto timedResult = _serviceManager->getSearchResult(searchId);
        if (!timedResult) {
            res.code = 404;
            res.body = "Query not found or expired";
            return res;
        }
        
        if (!timedResult->completed) {
            crow::json::wvalue response;
            response["status"] = "pending";
            response["searchId"] = searchId;
            res.code = 202; // Accepted
            res.body = response.dump();
            res.add_header("Content-Type", "application/json");
            return res;
        }
        
        if (!timedResult->error.empty()) {
            crow::json::wvalue response;
            response["status"] = "error";
            response["error"] = timedResult->error;
            res.code = 500;
            res.body = response.dump();
            res.add_header("Content-Type", "application/json");
            return res;
        }
        
        if (!timedResult->data) {
            res.code = 500;
            res.body = "No search data available";
            return res;
        }
        
        auto jsonResponse = resultToJson(*timedResult->data, offset, offset + limit);
        jsonResponse["searchId"] = searchId;
        jsonResponse["status"] = "completed";
        
        res.code = 200;
        res.body = jsonResponse.dump();
        res.add_header("Content-Type", "application/json");
        
    } catch (const std::exception& e) {
        LOG(ERROR) << "Error in handleGetQuery: " << e.what();
        res.code = 400;
        res.body = "Invalid parameter format: " + std::string(e.what());
    }
    
    return res;
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
crow::response WebService<SearchEngineT>::handleFastQuery(const crow::request& req) {
    crow::response res;
    res.add_header("Access-Control-Allow-Origin", "*");
    res.add_header("Content-Type", "application/json");
    
    try {
        auto json = crow::json::load(req.body);
        if (!json || !json.has("smiles") || !json.has("collectionId")) {
            crow::json::wvalue response;
            response["error"] = "Missing or invalid 'smiles' or 'collectionId' field";
            res.code = 400;
            res.body = response.dump();
            return res;
        }
        
        std::string smiles = json["smiles"].s();
        std::string collectionId = json["collectionId"].s();
        auto limit = json.has("limit") ? static_cast<size_t>(json["limit"].i()) : 10;
        
        if (!_serviceManager->isCollectionReady(collectionId)) {
            crow::json::wvalue response;
            response["error"] = "Collection not ready: " + collectionId;
            response["code"] = -1;
            res.code = 503;
            res.body = response.dump();
            return res;
        }
        
        auto result = _serviceManager->quickSearch(collectionId, smiles, limit);
        if (!result) {
            crow::json::wvalue response;
            response["error"] = "Search failed";
            response["code"] = -1;
            res.code = 500;
            res.body = response.dump();
            return res;
        }
        
        auto jsonResponse = resultToJson(*result, 0, limit);
        res.code = 200;
        res.body = jsonResponse.dump();
        
        LOG(INFO) << "FastQuery completed for SMILES: " << smiles << " in collection: " << collectionId << " with " << result->size() << " results";
        
    } catch (const std::exception& e) {
        LOG(ERROR) << "Error in handleFastQuery: " << e.what();
        crow::json::wvalue response;
        response["error"] = "Internal server error: " + std::string(e.what());
        response["code"] = -1;
        res.code = 500;
        res.body = response.dump();
    }
    
    return res;
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
crow::response WebService<SearchEngineT>::handleStatus(const crow::request& req) {
    crow::response res;
    res.add_header("Access-Control-Allow-Origin", "*");
    res.add_header("Content-Type", "application/json");
    
    crow::json::wvalue response;
    response["status"] = "running";
    auto collections = _serviceManager->getCollectionIds();
    response["collections_count"] = collections.size();
    
    res.code = 200;
    res.body = response.dump();
    
    return res;
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
crow::response WebService<SearchEngineT>::handleSmilesById(const crow::request& req) {
    crow::response res;
    res.add_header("Access-Control-Allow-Origin", "*");
    res.add_header("Content-Type", "application/json");
    
    try {
        const char* idStr = req.url_params.get("id");
        const char* collectionIdStr = req.url_params.get("collectionId");
        
        if (!idStr || !collectionIdStr) {
            crow::json::wvalue response;
            response["error"] = "Missing 'id' or 'collectionId' parameter";
            res.code = 400;
            res.body = response.dump();
            return res;
        }
        
        std::string collectionId = collectionIdStr;
        if (!_serviceManager->isCollectionReady(collectionId)) {
            crow::json::wvalue response;
            response["error"] = "Collection not ready: " + collectionId;
            res.code = 503;
            res.body = response.dump();
            return res;
        }
        
        size_t moleculeId = std::stoull(idStr);
        
        // Get SMILES string from the collection
        std::string smiles = _serviceManager->getSmilesById(collectionId, moleculeId);
        
        if (smiles.empty()) {
            crow::json::wvalue response;
            response["error"] = "Molecule ID not found in collection";
            res.code = 404;
            res.body = response.dump();
            return res;
        }
        
        crow::json::wvalue response;
        response["id"] = moleculeId;
        response["smiles"] = smiles;
        response["collectionId"] = collectionId;
        
        res.code = 200;
        res.body = response.dump();
        
        LOG(INFO) << "Retrieved SMILES for molecule ID " << moleculeId << " in collection " << collectionId << ": " << smiles;
        
    } catch (const std::invalid_argument& e) {
        LOG(ERROR) << "Invalid molecule ID format: " << e.what();
        crow::json::wvalue response;
        response["error"] = "Invalid ID format";
        res.code = 400;
        res.body = response.dump();
    } catch (const std::exception& e) {
        LOG(ERROR) << "Error in handleSmilesById: " << e.what();
        crow::json::wvalue response;
        response["error"] = "Internal server error: " + std::string(e.what());
        res.code = 500;
        res.body = response.dump();
    }
    
    return res;
}

// Collection management handlers
template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
crow::response WebService<SearchEngineT>::handleCreateCollection(const crow::request& req) {
    crow::response res;
    res.add_header("Access-Control-Allow-Origin", "*");
    res.add_header("Content-Type", "application/json");
    
    try {
        auto collectionId = _serviceManager->createCollection();
        
        crow::json::wvalue response;
        response["collectionId"] = collectionId;
        response["status"] = "created";
        
        res.code = 201;
        res.body = response.dump();
        
        LOG(INFO) << "Created collection: " << collectionId;
        
    } catch (const std::exception& e) {
        LOG(ERROR) << "Error in handleCreateCollection: " << e.what();
        crow::json::wvalue response;
        response["error"] = "Internal server error: " + std::string(e.what());
        res.code = 500;
        res.body = response.dump();
    }
    
    return res;
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
crow::response WebService<SearchEngineT>::handleAddMolecule(const crow::request& req, const std::string& collectionId) {
    crow::response res;
    res.add_header("Access-Control-Allow-Origin", "*");
    res.add_header("Content-Type", "application/json");
    
    try {
        auto json = crow::json::load(req.body);
        if (!json || !json.has("smiles")) {
            crow::json::wvalue response;
            response["error"] = "Missing 'smiles' field";
            res.code = 400;
            res.body = response.dump();
            return res;
        }
        
        std::string smiles = json["smiles"].s();
        
        if (!_serviceManager->collectionExists(collectionId)) {
            crow::json::wvalue response;
            response["error"] = "Collection not found: " + collectionId;
            res.code = 404;
            res.body = response.dump();
            return res;
        }
        
        bool success = _serviceManager->addMolecule(collectionId, smiles);
        if (!success) {
            crow::json::wvalue response;
            response["error"] = "Failed to add molecule. Collection may not be in building state.";
            res.code = 400;
            res.body = response.dump();
            return res;
        }
        
        crow::json::wvalue response;
        response["status"] = "added";
        response["collectionId"] = collectionId;
        response["collectionSize"] = _serviceManager->getCollectionSize(collectionId);
        
        res.code = 200;
        res.body = response.dump();
        
        LOG(INFO) << "Added molecule to collection " << collectionId;
        
    } catch (const std::exception& e) {
        LOG(ERROR) << "Error in handleAddMolecule: " << e.what();
        crow::json::wvalue response;
        response["error"] = "Internal server error: " + std::string(e.what());
        res.code = 500;
        res.body = response.dump();
    }
    
    return res;
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
crow::response WebService<SearchEngineT>::handleUploadCSV(const crow::request& req, const std::string& collectionId) {
    crow::response res;
    res.add_header("Access-Control-Allow-Origin", "*");
    res.add_header("Content-Type", "application/json");
    
    try {
        if (!_serviceManager->collectionExists(collectionId)) {
            crow::json::wvalue response;
            response["error"] = "Collection not found: " + collectionId;
            res.code = 404;
            res.body = response.dump();
            return res;
        }
        
        // Extract CSV content from multipart form data or raw body
        std::string csvContent;
        
        // Check if it's multipart form data
        auto contentType = req.get_header_value("Content-Type");
        if (contentType.find("multipart/form-data") != std::string::npos) {
            // Parse multipart form data to extract file content
            size_t boundaryPos = contentType.find("boundary=");
            if (boundaryPos == std::string::npos) {
                crow::json::wvalue response;
                response["error"] = "Invalid multipart form data - no boundary found";
                res.code = 400;
                res.body = response.dump();
                return res;
            }
            
            std::string boundary = "--" + contentType.substr(boundaryPos + 9);
            std::string body = req.body;
            
            // Find the file content between boundaries
            size_t startPos = body.find(boundary);
            if (startPos != std::string::npos) {
                startPos = body.find("\r\n\r\n", startPos);
                if (startPos != std::string::npos) {
                    startPos += 4; // Skip "\r\n\r\n"
                    size_t endPos = body.find("\r\n" + boundary, startPos);
                    if (endPos != std::string::npos) {
                        csvContent = body.substr(startPos, endPos - startPos);
                    }
                }
            }
        } else {
            // Assume raw CSV content in body
            csvContent = req.body;
        }
        
        if (csvContent.empty()) {
            crow::json::wvalue response;
            response["error"] = "No CSV content found in request";
            res.code = 400;
            res.body = response.dump();
            return res;
        }
        
        // Parse CSV content directly from string
        std::vector<std::string> smilesList;
        std::istringstream csvStream(csvContent);
        std::string line;
        
        while (std::getline(csvStream, line)) {
            if (line.empty()) continue;
            
            std::istringstream lineStream(line);
            std::string id, smiles;
            lineStream >> id >> smiles;
            
            if (!smiles.empty()) {
                smilesList.push_back(smiles);
            } else {
                LOG(WARNING) << "Skipped malformed CSV line: " << line;
            }
        }
        
        if (smilesList.empty()) {
            crow::json::wvalue response;
            response["error"] = "No valid SMILES found in CSV";
            res.code = 400;
            res.body = response.dump();
            return res;
        }
        
        // Add molecules to the collection using SmilesStorage
        bool success = _serviceManager->addMolecules(collectionId, smilesList);
        if (!success) {
            crow::json::wvalue response;
            response["error"] = "Failed to add molecules. Collection may not be in building state.";
            res.code = 400;
            res.body = response.dump();
            return res;
        }
        
        crow::json::wvalue response;
        response["status"] = "uploaded";
        response["collectionId"] = collectionId;
        response["addedCount"] = smilesList.size();
        response["collectionSize"] = _serviceManager->getCollectionSize(collectionId);
        
        res.code = 200;
        res.body = response.dump();
        
        LOG(INFO) << "Uploaded " << smilesList.size() << " molecules from CSV to collection " << collectionId;
        
    } catch (const std::exception& e) {
        LOG(ERROR) << "Error in handleUploadCSV: " << e.what();
        crow::json::wvalue response;
        response["error"] = "Internal server error: " + std::string(e.what());
        res.code = 500;
        res.body = response.dump();
    }
    
    return res;
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
crow::response WebService<SearchEngineT>::handleBuildIndex(const crow::request& req, const std::string& collectionId) {
    crow::response res;
    res.add_header("Access-Control-Allow-Origin", "*");
    res.add_header("Content-Type", "application/json");
    
    try {
        
        if (!_serviceManager->collectionExists(collectionId)) {
            crow::json::wvalue response;
            response["error"] = "Collection not found: " + collectionId;
            res.code = 404;
            res.body = response.dump();
            return res;
        }
        
        bool success = _serviceManager->buildIndex(collectionId);
        if (!success) {
            auto error = _serviceManager->getCollectionError(collectionId);
            crow::json::wvalue response;
            response["error"] = "Failed to start index building: " + error;
            res.code = 400;
            res.body = response.dump();
            return res;
        }
        
        crow::json::wvalue response;
        response["status"] = "building";
        response["collectionId"] = collectionId;
        response["message"] = "Index building started";
        
        res.code = 202; // Accepted
        res.body = response.dump();
        
        LOG(INFO) << "Started building index for collection " << collectionId;
        
    } catch (const std::exception& e) {
        LOG(ERROR) << "Error in handleBuildIndex: " << e.what();
        crow::json::wvalue response;
        response["error"] = "Internal server error: " + std::string(e.what());
        res.code = 500;
        res.body = response.dump();
    }
    
    return res;
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
crow::response WebService<SearchEngineT>::handleGetCollections(const crow::request& req) {
    crow::response res;
    res.add_header("Access-Control-Allow-Origin", "*");
    res.add_header("Content-Type", "application/json");
    
    try {
        auto collectionIds = _serviceManager->getCollectionIds();
        
        crow::json::wvalue response;
        crow::json::wvalue collectionsArray = crow::json::wvalue::list();
        
        for (size_t i = 0; i < collectionIds.size(); ++i) {
            const auto& id = collectionIds[i];
            crow::json::wvalue collection;
            collection["id"] = id;
            collection["status"] = static_cast<int>(_serviceManager->getCollectionStatus(id));
            collection["size"] = _serviceManager->getCollectionSize(id);
            collectionsArray[i] = std::move(collection);
        }
        
        response["collections"] = std::move(collectionsArray);
        response["count"] = collectionIds.size();
        
        res.code = 200;
        res.body = response.dump();
        
    } catch (const std::exception& e) {
        LOG(ERROR) << "Error in handleGetCollections: " << e.what();
        crow::json::wvalue response;
        response["error"] = "Internal server error: " + std::string(e.what());
        res.code = 500;
        res.body = response.dump();
    }
    
    return res;
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
crow::response WebService<SearchEngineT>::handleGetCollectionInfo(const crow::request& req, const std::string& collectionId) {
    crow::response res;
    res.add_header("Access-Control-Allow-Origin", "*");
    res.add_header("Content-Type", "application/json");
    
    try {
        
        if (!_serviceManager->collectionExists(collectionId)) {
            crow::json::wvalue response;
            response["error"] = "Collection not found: " + collectionId;
            res.code = 404;
            res.body = response.dump();
            return res;
        }
        
        crow::json::wvalue response;
        response["id"] = collectionId;
        response["status"] = static_cast<int>(_serviceManager->getCollectionStatus(collectionId));
        response["size"] = _serviceManager->getCollectionSize(collectionId);
        
        auto error = _serviceManager->getCollectionError(collectionId);
        if (!error.empty()) {
            response["error"] = error;
        }
        
        res.code = 200;
        res.body = response.dump();
        
    } catch (const std::exception& e) {
        LOG(ERROR) << "Error in handleGetCollectionInfo: " << e.what();
        crow::json::wvalue response;
        response["error"] = "Internal server error: " + std::string(e.what());
        res.code = 500;
        res.body = response.dump();
    }
    
    return res;
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
crow::json::wvalue WebService<SearchEngineT>::resultToJson(const SearchResult<size_t>& result, size_t offset, size_t limit) {
    crow::json::wvalue response;
    
    response["total"] = result.size();
    response["offset"] = offset;
    
    crow::json::wvalue results = crow::json::wvalue::list();
    
    size_t end = std::min(limit, result.size());
    size_t actualOffset = std::min(offset, result.size());
    
    for (size_t i = actualOffset; i < end && i < actualOffset + (limit - offset); ++i) {
        crow::json::wvalue item;
        item["smiles"] = result.get(i);
        item["index"] = i;
        results[i - actualOffset] = std::move(item);
    }
    
    response["results"] = std::move(results);
    response["count"] = end - actualOffset;
    
    return response;
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
std::pair<size_t, size_t> WebService<SearchEngineT>::extractBounds(const crow::json::rvalue& json) {
    size_t offset = json.has("offset") ? static_cast<size_t>(json["offset"].i()) : 0;
    size_t limit = json.has("limit") ? static_cast<size_t>(json["limit"].i()) : 100;
    return {offset, limit};
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
void WebService<SearchEngineT>::startCleanupThread() {
    _cleanupThread = std::make_unique<std::thread>([this]() {
        while (!_shutdownFlag) {
            std::this_thread::sleep_for(std::chrono::seconds(60));
            
            if (!_shutdownFlag) {
                try {
                    _serviceManager->cleanupExpiredQueries();
                } catch (const std::exception& e) {
                    LOG(ERROR) << "Error during cleanup: " << e.what();
                }
            }
        }
    });
} 