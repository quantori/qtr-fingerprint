#pragma once

#include <utility>

#include "RunMode.h"
#include "IndigoRamFilter.h"
#include "IndigoDriveFilter.h"

namespace qtr {
    class FromFileMode : public RunMode {
    private:
        std::filesystem::path _inputFile;

    public:
        inline FromFileMode(std::shared_ptr<SearchData> searchData, std::filesystem::path inputFile)
                : RunMode(std::move(searchData)), _inputFile(std::move(inputFile)) {}


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
                this->_searchData->timeTicker.tick();
                auto queryData = this->_searchData->search(queries[i], PropertiesFilter::Bounds());
                if (queryData == nullptr) {
                    ++skipped;
                    continue;
                }
                queryData->waitAllTasks();
                LOG(INFO) << "Found " << queryData->getCurrentAnswersCount() << " answers";
                times.emplace_back(
                        this->_searchData->timeTicker.tick("search molecule " + std::to_string(i) + ": " + queries[i]));
            }


            double mean = std::accumulate(times.begin(), times.end(), 0.0) / double(times.size());
            double min = *std::min_element(times.begin(), times.end());
            double max = *std::max_element(times.begin(), times.end());
            std::sort(times.begin(), times.end());
            double median = times[times.size() / 2];
            static const int percentiles[] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 95};
            std::ostringstream percentilesStatStream;
            for (int i = 0; i < std::size(percentiles); i++) {
                percentilesStatStream << percentiles[i] << "%: " << times[times.size() * percentiles[i] / 100];
                if (i + 1 != std::size(percentiles))
                    percentilesStatStream << " | ";
            }

            LOG(INFO) << "skipped queries: " << skipped;
            LOG(INFO) << "  mean: " << mean;
            LOG(INFO) << "   max: " << max;
            LOG(INFO) << "   min: " << min;
            LOG(INFO) << "median: " << median;
            LOG(INFO) << percentilesStatStream.str();
            LOG(INFO) << "Total search time: " << BallTreeSearchEngine::ballTreeSearchTimer;
            LOG(INFO) << "Total indigo time: " << IndigoFilter::indigoFilteringTimer;
            LOG(INFO) << "indigo percentage: "
                      << IndigoFilter::indigoFilteringTimer / BallTreeSearchEngine::ballTreeSearchTimer * 100 << "%";
            LOG(INFO) << "overdue queries: " << BallTreeQueryData::timedOutCounter;
        }

    };
} // qtr