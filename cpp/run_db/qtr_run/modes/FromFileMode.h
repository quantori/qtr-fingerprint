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
        bool _fingerprintProvided;

    public:
        inline FromFileMode(std::shared_ptr<SearchData> searchData, std::filesystem::path inputFile,
                            std::filesystem::path summaryFile, bool fingerprintProvided)
                : RunMode(std::move(searchData)), _inputFile(std::move(inputFile)),
                  _summaryFile(std::move(summaryFile)), _fingerprintProvided(fingerprintProvided) {}

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

            out << "skipped queries: " << skippedQueries << '\n';
            out << "  mean: " << mean << '\n';
            out << "   max: " << max << '\n';
            out << "   min: " << min << '\n';
            out << "median: " << median << '\n';
            out << percentilesStatStream.str() << '\n';
            out << "overdue queries: " << BallTreeQueryData::timedOutCounter << '\n';
        }

        inline static std::vector<SearchData::Query>
        loadQueriesFromFile(const std::filesystem::path &inputFile, bool fingerprintProvided) {
            std::ifstream input(inputFile);
            std::vector<SearchData::Query> queries;
            while (input.peek() != EOF) {
                auto smiles = std::make_unique<std::string>();
                input >> *smiles;
                std::unique_ptr<Fingerprint> fingerprint = nullptr;
                if (fingerprintProvided) {
                    std::string fingerprintStr;
                    input >> fingerprintStr;
                    fingerprint = std::make_unique<Fingerprint>(fingerprintStr);
                }
                std::string otherInfoInLine;
                std::getline(input, otherInfoInLine);
                queries.emplace_back(std::move(smiles), std::move(fingerprint));
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
            auto queries = loadQueriesFromFile(_inputFile, _fingerprintProvided);
            std::vector<float> times;
            std::vector<size_t> answerCounters;
            size_t skipped = 0;
            for (size_t i = 0; i < queries.size(); i++) {
                assert(queries[i].smiles != nullptr);
                LOG(INFO) << "Start search for " << i << ": " << *queries[i].smiles;
                ProfilingTimer profilingTimer("Query processing");
                auto queryData = this->_searchData->search(queries[i], PropertiesFilter::Bounds());
                if (queryData == nullptr) {
                    ++skipped;
                    continue;
                }
                queryData->waitAllTasks();
                auto ansCount = queryData->getCurrentAnswersCount();
                auto answers = queryData->getAnswers(0, std::min((size_t) 10, ansCount)).second;
                answerCounters.emplace_back(ansCount);
                LOG(INFO) << "Found " << ansCount << " answers. Examples: " << [&answers]() {
                    std::stringstream s;
                    for (auto &ans: answers)
                        s << ans << ' ';
                    return s.str();
                }();
                auto queryDuration = profilingTimer.stop();
                LOG(INFO) << queryDuration << " seconds spent to process molecule " << i << ": " << *queries[i].smiles;
                times.emplace_back(queryDuration);
            }

            showStatistics(times, skipped, std::cout);
            showSummary(answerCounters, _summaryFile);
            ProfilingPool::showStatistics(std::cout);
        }


    };
} // qtr