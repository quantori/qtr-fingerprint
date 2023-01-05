#include "WebMode.h"

#include <utility>

using namespace std;

namespace qtr {


    WebMode::WebMode(const BallTreeSearchEngine &ballTree, shared_ptr<const SmilesTable> smilesTable,
                     TimeTicker &timeTicker,
                     uint64_t ansCount, uint64_t threadsCount,
                     std::filesystem::path &idToStringDirPath) : _ballTree(ballTree),
                                                                 _smilesTable(std::move(smilesTable)),
                                                                 _ansCount(ansCount),
                                                                 _threadsCount(threadsCount),
                                                                 _idConverter(idToStringDirPath) {}


    crow::json::wvalue
    WebMode::prepareResponse(BallTreeQueryData &queryData, size_t minOffset, size_t maxOffset) {
        crow::json::wvalue::list response;
        cout << minOffset << " " << maxOffset << endl;
        auto [isFinish, result] = queryData.getAnswers(minOffset, maxOffset);
        LOG(INFO) << "found " << result.size() << " results";
        response.reserve(result.size());
        response.emplace_back(crow::json::wvalue{{"isFinished", isFinish}});
        for (auto res: result) {
            auto [id, libraryId] = _idConverter.fromDbId(res);
            response.emplace_back(crow::json::wvalue{{"id",        id},
                                                     {"libraryId", libraryId}});
        }

        return crow::json::wvalue{response};
    }


    void WebMode::run() {
        crow::SimpleApp app;
        mutex newTaskMutex;

        vector<string> smiles;
        vector<unique_ptr<BallTreeQueryData>> queries;
        unordered_map<string, uint64_t> queryToId;

        CROW_ROUTE(app, "/query").methods(crow::HTTPMethod::POST)(
                [&queryToId, &queries, &newTaskMutex, this, &smiles](
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
                        auto [error, queryData] = doSearch(currSmiles, _ballTree, _smilesTable, _ansCount,
                                                           _threadsCount);
                        queries.emplace_back(std::move(queryData));
                    }
                    return crow::response(to_string(queryToId[currSmiles]));
                });

        CROW_ROUTE(app, "/query").methods(crow::HTTPMethod::GET)(
                [&queries, this](const crow::request &req) {
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