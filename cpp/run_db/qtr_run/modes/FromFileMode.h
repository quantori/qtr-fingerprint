#pragma once

#include <utility>

#include "RunMode.h"
#include "IndigoSmilesRamFilter.h"
#include "IndigoSmilesDriveFilter.h"
#include "Profiling.h"

namespace qtr {
    class FromFileMode : public RunMode {
    private:
        std::filesystem::path _inputFile;
        std::filesystem::path _summaryFile;

    public:
        inline FromFileMode(std::shared_ptr<SearchData> searchData, std::filesystem::path inputFile,
                            std::filesystem::path summaryFile)
                : RunMode(std::move(searchData)), _inputFile(std::move(inputFile)), _summaryFile(summaryFile) {}

        inline static void showStatistics(std::vector<float> times, size_t skippedQueries, std::ostream &out) {
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

            out << "skipped queries: " << skippedQueries;
            out << "  mean: " << mean;
            out << "   max: " << max;
            out << "   min: " << min;
            out << "median: " << median;
            out << percentilesStatStream.str();
            out << "Total search time: " << BallTreeSearchEngine::ballTreeSearchTimer;
            out << "Total indigo time: " << IndigoSmilesFilter::indigoFilteringTimer;
            out << "indigo percentage: "
                << IndigoSmilesFilter::indigoFilteringTimer / BallTreeSearchEngine::ballTreeSearchTimer * 100
                << "%";
            out << "overdue queries: " << BallTreeQueryData::timedOutCounter;
        }

        inline static std::vector<std::string> loadQueriesFromFile(const std::filesystem::path &inputFile) {
            std::ifstream input(inputFile);
            std::vector<std::string> queries;
            while (input.peek() != EOF) {
                std::string query;
                input >> query;
                std::string otherInfoInLine;
                std::getline(input, otherInfoInLine);
                queries.emplace_back(query);
            }
            LOG(INFO) << "Loaded " << queries.size() << " queries";
            return queries;
        }

        inline static void
        showSummary(const std::vector<size_t> &answerCounters, const std::filesystem::path &summaryFile) {
            if (summaryFile.empty())
                return;
            std::ofstream out(summaryFile);
            for (float i: answerCounters) {
                out << i << '\n';
            }
        }

        inline void run() override {
            auto queries = loadQueriesFromFile(_inputFile);
            std::vector<float> times;
            std::vector<size_t> answerCounters;
            size_t skipped = 0;
            for (size_t i = 0; i < queries.size(); i++) {
                LOG(INFO) << "Start search for " << i << ": " << queries[i];
                ProfilingTimer profilingTimer("Query processing");
                auto queryData = this->_searchData->search(queries[i], PropertiesFilter::Bounds());
                if (queryData == nullptr) {
                    ++skipped;
                    continue;
                }
                queryData->waitAllTasks();
                auto ansCount = queryData->getCurrentAnswersCount();
                answerCounters.emplace_back(ansCount);
                LOG(INFO) << "Found " << ansCount << " answers";
                auto queryDuration = profilingTimer.stop();
                LOG(INFO) << queryDuration << " seconds spent to process molecule " << i << ": " << queries[i];
                times.emplace_back(queryDuration);
            }

            showStatistics(times, skipped, std::cout);
            showSummary(answerCounters, _summaryFile);
        }


    };
} // qtr