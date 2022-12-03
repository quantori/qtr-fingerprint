#include "WebMode.h"

using namespace qtr;

namespace qtr {


    WebMode::WebMode(const BallTreeSearchEngine &ballTree, const SmilesTable &smilesTable, TimeTicker &timeTicker,
                     uint64_t ansCount, uint64_t startSearchDepth,
                     std::filesystem::path &idToStringDirPath) : ballTree(ballTree),
                                                                 smilesTable(smilesTable),
                                                                 ansCount(ansCount),
                                                                 startSearchDepth(startSearchDepth),
                                                                 idConverter(idToStringDirPath) {}


    crow::json::wvalue
    WebMode::prepareResponse(const std::vector<uint64_t> &ids, size_t minOffset, size_t maxOffset) {
        crow::json::wvalue::list response;
        std::cout << minOffset << " " << maxOffset << std::endl;
        if (!ids.empty()) {
            auto rightBorder = std::min(ids.size(), maxOffset);
            auto leftBorder = std::min(ids.size(), minOffset);
            response.reserve(rightBorder - leftBorder);
            for (size_t i = leftBorder; i < rightBorder; ++i) {
                auto [id, libraryId] = idConverter.fromDbId(ids[i]);
                response.emplace_back(crow::json::wvalue{{"id",        id},
                                                         {"libraryId", libraryId}});
            }
        }
        return crow::json::wvalue{response};
    }


    void WebMode::run() {
        crow::SimpleApp app;
        std::mutex newTaskMutex;
        uint64_t _queryIdTicker = 0;

        std::unordered_map<uint64_t, std::future<std::pair<bool, std::vector<uint64_t>>>> tasks;
        std::unordered_map<uint64_t, std::vector<uint64_t>> resultTable;
        std::unordered_map<std::string, uint64_t> queryToId;

        CROW_ROUTE(app, "/query").methods(crow::HTTPMethod::POST)(
                [&tasks, &_queryIdTicker, &queryToId, &newTaskMutex, this](
                        const crow::request &req) {
                    auto body = crow::json::load(req.body);
                    if (!body)
                        return crow::response(400);
                    std::cout << body["smiles"] << std::endl;
                    std::string smiles = body["smiles"].s();
                    std::lock_guard<std::mutex> lock(newTaskMutex);
                    if (!queryToId.contains(smiles)) {
                        _queryIdTicker += 1;
                        tasks[_queryIdTicker] = std::async(std::launch::async, doSearch, smiles,
                                                           std::ref(ballTree), std::ref(smilesTable), ansCount,
                                                           startSearchDepth);
                        queryToId[smiles] = _queryIdTicker;
                    }
                    return crow::response(std::to_string(queryToId[smiles]));
                });

        CROW_ROUTE(app, "/query").methods(crow::HTTPMethod::GET)(
                [&tasks, &resultTable, this](const crow::request &req) {
                    auto searchId = crow::utility::lexical_cast<uint64_t>(req.url_params.get("searchId"));
                    auto offset = crow::utility::lexical_cast<int>(req.url_params.get("offset"));
                    auto limit = crow::utility::lexical_cast<int>(req.url_params.get("limit"));
                    if (resultTable.contains(searchId))
                        return prepareResponse(resultTable[searchId], offset, offset + limit);
                    if (!tasks.contains(searchId))
                        return crow::json::wvalue();
                    const auto &task = tasks[searchId];
                    if (task.wait_for(std::chrono::seconds(0)) == std::future_status::ready) {
                        auto [isSkipped, result] = tasks[searchId].get();
                        resultTable[searchId] = std::move(result);
                    } else {
                        return crow::json::wvalue();
                    }
                    return prepareResponse(resultTable[searchId], offset, offset + limit);
                });

        app.port(8080).concurrency(2).run();

    }
} // namespace qtr