#include "WebMode.h"
#include "search_data/QtrSearchData.h"

#include <utility>

using namespace std;

namespace qtr {

    namespace {
        pair<string, string> convertId(CIDType id, const shared_ptr<SearchData>& searchData) {
            auto qtrSearchData = dynamic_cast<const QtrSearchData *>(searchData.get());
            if (qtrSearchData != nullptr)
                return qtrSearchData->idConverter->fromDbId(id);
            else
                return {to_string(id), ""};
        }
    }


    WebMode::WebMode(std::shared_ptr<SearchData> searchData) : RunMode(std::move(searchData)) {}

    static PropertiesFilter::Bounds extractBounds(const crow::json::rvalue &json) {
        PropertiesFilter::Bounds bounds;

        for (size_t i = 0; i < std::size(PropertiesFilter::propertyNames); ++i) {
            auto &propertyName = PropertiesFilter::propertyNames[i];
            if (!json.has(propertyName))
                continue;
            auto &jsonBounds = json[propertyName];

            if (jsonBounds.has("min")) {
                bounds.minBounds[i] = (float) jsonBounds["min"].d();
                LOG(INFO) << "Set min bound " << propertyName << ": " << bounds.minBounds[i];
            }
            if (jsonBounds.has("max")) {
                bounds.maxBounds[i] = (float) jsonBounds["max"].d();
                LOG(INFO) << "Set max bound " << propertyName << ": " << bounds.maxBounds[i];
            }
        }

        return bounds;
    }

    void WebMode::run() {
        crow::SimpleApp app;
        mutex newTaskMutex;

        vector<string> smilesList;
        vector<unique_ptr<QueryData<CIDType>>> queries;

        CROW_ROUTE(app, "/query")
                .methods(crow::HTTPMethod::POST)(
                        [&queries, &newTaskMutex, this, &smilesList](const crow::request &req) {
                            auto json = crow::json::load(req.body);
                            if (!json)
                                return crow::response(400);

                            lock_guard<mutex> lock(newTaskMutex);
                            smilesList.emplace_back(json["smiles"].s());
                            string &currSmiles = smilesList.back();

                            auto bounds = extractBounds(json);
                            size_t queryId = queries.size();
                            auto queryData = _searchData->search(currSmiles, bounds);
                            if (queryData == nullptr) {
                                LOG(WARNING) << "Cannot start search for smilesList: " << currSmiles;
                                return crow::response(to_string(-1));
                            }
                            queries.emplace_back(std::move(queryData));

                            return crow::response(to_string(queryId));
                        });

        CROW_ROUTE(app, "/query")
                .methods(crow::HTTPMethod::GET)([&queries, this](const crow::request &req) {
                    auto searchId = crow::utility::lexical_cast<uint64_t>(req.url_params.get("searchId"));
                    auto offset = crow::utility::lexical_cast<uint64_t>(req.url_params.get("offset"));
                    auto limit = crow::utility::lexical_cast<uint64_t>(req.url_params.get("limit"));

                    if (searchId >= queries.size())
                        return crow::json::wvalue();

                    return prepareResponse(*queries[searchId], offset, offset + limit);
                });

        CROW_ROUTE(app, "/fastquery")
                .methods(crow::HTTPMethod::GET)([&newTaskMutex, this](const crow::request &req) {
                    auto limit = crow::utility::lexical_cast<uint64_t>(req.url_params.get("limit"));
                    auto json = crow::json::load(req.body);
                    if (!json)
                        return crow::json::wvalue(to_string(-1));
                    auto smiles = json["smiles"].s();
                    auto bounds = extractBounds(json);

                    lock_guard<mutex> lock(newTaskMutex);
                    auto queryData = _searchData->search(smiles, bounds);
                    if (queryData == nullptr) {
                        LOG(WARNING) << "Cannot start search for smilesList: " << smiles;
                        return crow::json::wvalue(to_string(-1));
                    }
                    queryData->waitAllTasks();
                    return prepareResponse(*queryData, 0, limit);
                });

        app.port(8080).concurrency(2).run();

    }

    crow::json::wvalue WebMode::prepareResponse(QueryData<CIDType> &queryData, size_t minOffset, size_t maxOffset) {
        cout << minOffset << " " << maxOffset << endl;
        auto [isFinished, result] = queryData.getAnswers(minOffset, maxOffset);
        LOG(INFO) << "found " << result.size() << " results";
        crow::json::wvalue::list molecules;
        molecules.reserve(result.size());
        for (auto res: result) {
            auto [id, libraryId] = convertId(res, _searchData);
            molecules.emplace_back(crow::json::wvalue{{"id",        id},
                                                      {"libraryId", libraryId}});
        }
        return crow::json::wvalue{
                {"molecules",  molecules},
                {"isFinished", isFinished}
        };
    }

} // namespace qtr