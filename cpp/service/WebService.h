#pragma once

#include <memory>
#include <atomic>
#include <thread>
#include <chrono>
#include <sstream>

#include <crow.h>
#include <glog/logging.h>

#include "ServiceManager.h"
#include "utils/Config.h"

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
}

template<typename SearchEngineT>
requires SearchEngineInterface<SearchEngineT>
crow::response WebService<SearchEngineT>::handlePostQuery(const crow::request& req) {
    crow::response res;
    res.add_header("Access-Control-Allow-Origin", "*");
    
    try {
        auto json = crow::json::load(req.body);
        if (!json || !json.has("smiles")) {
            res.code = 400;
            res.body = "Missing or invalid 'smiles' field";
            return res;
        }
        
        auto smiles = json["smiles"].s();
        auto maxResults = json.has("maxResults") ? static_cast<size_t>(json["maxResults"].i()) : 1000;
        auto timeLimit = json.has("timeLimit") ? json["timeLimit"].d() : 60.0;
        
        if (!_serviceManager->isEngineReady()) {
            res.code = 503;
            res.body = "Search engine not ready";
            return res;
        }
        
        auto queryId = _serviceManager->submitSearch(smiles, maxResults, timeLimit);
        
        res.code = 200;
        res.body = std::to_string(queryId);
        
        LOG(INFO) << "Submitted query with ID: " << queryId << " for SMILES: " << smiles;
        
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
        if (!json || !json.has("smiles")) {
            crow::json::wvalue response;
            response["error"] = "Missing or invalid 'smiles' field";
            res.code = 400;
            res.body = response.dump();
            return res;
        }
        
        auto smiles = json["smiles"].s();
        auto limit = json.has("limit") ? static_cast<size_t>(json["limit"].i()) : 10;
        
        if (!_serviceManager->isEngineReady()) {
            crow::json::wvalue response;
            response["error"] = "Search engine not ready";
            response["code"] = -1;
            res.code = 503;
            res.body = response.dump();
            return res;
        }
        
        auto result = _serviceManager->quickSearch(smiles, limit);
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
        
        LOG(INFO) << "FastQuery completed for SMILES: " << smiles << " with " << result->size() << " results";
        
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
    response["engine_ready"] = _serviceManager->isEngineReady();
    
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
        
        if (!idStr) {
            crow::json::wvalue response;
            response["error"] = "Missing 'id' parameter";
            res.code = 400;
            res.body = response.dump();
            return res;
        }
        
        if (!_serviceManager->isEngineReady()) {
            crow::json::wvalue response;
            response["error"] = "Search engine not ready";
            res.code = 503;
            res.body = response.dump();
            return res;
        }
        
        size_t moleculeId = std::stoull(idStr);
        
        // Get SMILES string using the search engine's resultToSmiles method
        std::string smiles = _serviceManager->getSmilesById(moleculeId);
        
        if (smiles.empty()) {
            crow::json::wvalue response;
            response["error"] = "Molecule ID not found";
            res.code = 404;
            res.body = response.dump();
            return res;
        }
        
        crow::json::wvalue response;
        response["id"] = moleculeId;
        response["smiles"] = smiles;
        
        res.code = 200;
        res.body = response.dump();
        
        LOG(INFO) << "Retrieved SMILES for molecule ID " << moleculeId << ": " << smiles;
        
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