#include "WebMode.h"

using namespace qtr;
using namespace std;

namespace qtr {


    WebMode::WebMode(const BallTreeSearchEngine &ballTree, const SmilesTable &smilesTable, TimeTicker &timeTicker,
                     uint64_t ansCount, uint64_t startSearchDepth) : ballTree(ballTree),
                                                                     smilesTable(smilesTable),
                                                                     ansCount(ansCount),
                                                                     startSearchDepth(startSearchDepth) {}


    crow::json::wvalue
    WebMode::prepareResponse(BallTreeSearchEngine::QueryData &queryData, size_t minOffset, size_t maxOffset) {
        crow::json::wvalue::list response;
        cout << minOffset << " " << maxOffset << endl;
        lock_guard<mutex> lock(queryData.resultLock);
        LOG(INFO) << "found " << queryData.result.size() << " results";
        if (!queryData.result.empty()) {
            copy(queryData.result.begin() + min(queryData.result.size(), minOffset),
                 queryData.result.begin() + min(queryData.result.size(), maxOffset),
                 back_inserter(response));
        }
        return crow::json::wvalue{response};
    }


    void WebMode::run() {
        crow::SimpleApp app;
        mutex newTaskMutex;

        vector<vector<future<void>>> tasks;
        vector<string> smiles;
        vector<unique_ptr<BallTreeSearchEngine::QueryData>> queries;
        unordered_map<string, uint64_t> queryToId;

        CROW_ROUTE(app, "/query").methods(crow::HTTPMethod::POST)(
                [&queryToId, &queries, &tasks, &newTaskMutex, this, &smiles](
                        const crow::request &req) {
                    auto body = crow::json::load(req.body);
                    if (!body)
                        return crow::response(400);
                    cout << body["smiles"] << endl;
                    lock_guard<mutex> lock(newTaskMutex);
                    smiles.emplace_back(body["smiles"].s());
                    string &currSmiles = smiles.back();
                    if (!queryToId.contains(currSmiles)) {
                        queryToId[currSmiles] = queries.size();
                        queries.emplace_back(make_unique<BallTreeSearchEngine::QueryData>(ansCount));
                        tasks.resize(queries.size());
                        tasks.back() = doSearch(currSmiles,
                                                *queries.back(), ballTree,
                                                smilesTable, startSearchDepth).second;
                    }
                    return crow::response(to_string(queryToId[currSmiles]));
                });

        CROW_ROUTE(app, "/query").methods(crow::HTTPMethod::GET)(
                [&queries](const crow::request &req) {
                    auto searchId = crow::utility::lexical_cast<uint64_t>(req.url_params.get("searchId"));
                    auto offset = crow::utility::lexical_cast<int>(req.url_params.get("offset"));
                    auto limit = crow::utility::lexical_cast<int>(req.url_params.get("limit"));
                    if (searchId >= queries.size())
                        return crow::json::wvalue();
                    return prepareResponse(*queries[searchId], offset, offset + limit);
                });

        app.port(8080).concurrency(2).run();

    }
} // namespace qtr