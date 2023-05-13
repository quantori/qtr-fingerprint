#pragma once

#include <utility>

#include "RunMode.h"
#include "RunDbUtils.h"
#include "IndigoRamFilter.h"
#include "IndigoDriveFilter.h"

namespace qtr {
    class FromFileMode : public RunMode {
    private:
        std::shared_ptr<const SearchData> _searchData;
        std::filesystem::path _inputFile;

    public:
        inline FromFileMode(std::shared_ptr<const SearchData> searchData, std::filesystem::path inputFile)
                : _searchData(std::move(searchData)), _inputFile(std::move(inputFile)) {}


        inline void run() override {
            std::ifstream input(_inputFile);
            std::vector<std::string> queries;
            while (input.peek() != EOF) {
                std::string query;
                input >> query;
                std::string otherInfoInLine;
                std::getline(input, otherInfoInLine);
                queries.emplace_back(query);
            }
            std::vector<double> times;
            size_t skipped = 0;
            LOG(INFO) << "Loaded " << queries.size() << " queries";
            for (size_t i = 0; i < queries.size(); i++) {
                LOG(INFO) << "Start search for " << i << ": " << queries[i];
                _searchData->timeTicker.tick();
                auto [error, queryData] = runSearch(*_searchData, queries[i], PropertiesFilter::Bounds());
                if (error) {
                    ++skipped;
                    continue;
                }
                queryData->waitAllTasks();
                LOG(INFO) << "Found " << queryData->getCurrentAnswersCount() << " answers";
                times.emplace_back(
                        _searchData->timeTicker.tick("search molecule " + std::to_string(i) + ": " + queries[i]));
            }


            double mean = std::accumulate(times.begin(), times.end(), 0.0) / double(times.size());
            double min = *std::min_element(times.begin(), times.end());
            double max = *std::max_element(times.begin(), times.end());
            std::sort(times.begin(), times.end());
            double median = times[times.size() / 2];
            double p60 = times[times.size() * 6 / 10];
            double p70 = times[times.size() * 7 / 10];
            double p80 = times[times.size() * 8 / 10];
            double p90 = times[times.size() * 9 / 10];
            double p95 = times[times.size() * 95 / 100];
            LOG(INFO) << "skipped queries: " << skipped;
            LOG(INFO) << "  mean: " << mean;
            LOG(INFO) << "   max: " << max;
            LOG(INFO) << "   min: " << min;
            LOG(INFO) << "median: " << median;
            LOG(INFO) << "60%: " << p60 << " | 70%: " << p70 << " | 80%: " << p80 << " | 90%: " << p90 << " | 95%: "
                      << p95;
            LOG(INFO) << "Total search time: " << BallTreeSearchEngine::ballTreeSearchTimer;
            LOG(INFO) << "Total indigo time: "
                      << IndigoFilter::indigoFilteringTimer;
            LOG(INFO) << "indigo percentage: "
                      << IndigoFilter::indigoFilteringTimer /
                         BallTreeSearchEngine::ballTreeSearchTimer * 100 << "%";
            LOG(INFO) << "overdue queries: " << BallTreeQueryData::timedOutCounter;
        }

    };
} // qtr