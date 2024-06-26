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
        uint64_t _workers;

    public:
        inline FromFileMode(std::shared_ptr<SearchData> searchData, std::filesystem::path inputFile,
                            std::filesystem::path summaryFile, bool fingerprintProvided, uint64_t workers)
                : RunMode(std::move(searchData)), _inputFile(std::move(inputFile)),
                  _summaryFile(std::move(summaryFile)), _fingerprintProvided(fingerprintProvided), _workers(workers) {}

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
            out << "overdue queries: " << QueryDataWithFingerprint::timedOutCounter << '\n';
        }

        inline static std::vector<SearchData::Query>
        loadQueriesFromFile(const std::filesystem::path &inputFile, bool fingerprintProvided, BaseLibrary baseLibrary) {
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
                queries.emplace_back(std::move(smiles), std::move(fingerprint), baseLibrary);
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

        inline static std::vector<std::vector<SearchData::Query>>
        splitQueries(std::vector<SearchData::Query> &&queries, size_t count) {
            std::vector<std::vector<SearchData::Query>> result(count);
            for (size_t i = 0; i < queries.size(); i++) {
                result[i % count].push_back(std::move(queries[i]));
            }
            return result;
        }

        inline static void
        processGroup(const std::vector<SearchData::Query> &group, const std::shared_ptr<SearchData> &searchData,
                     size_t &skipped, std::mutex &mutex, size_t worker, size_t workers,
                     std::vector<size_t> &answerCounters, std::vector<float> &times) {
            for (size_t i = 0; i < group.size(); i++) {
                assert(group[i].smiles != nullptr);
                LOG(INFO) << "Start search for " << i * workers + worker << ": " << *group[i].smiles;
                ProfilingTimer profilingTimer("Query processing");

                auto queryData = searchData->search(group[i], PropertiesFilter::Bounds());
                if (queryData == nullptr) {
                    std::lock_guard<std::mutex> guard(mutex);
                    ++skipped;
                    continue;
                }

                queryData->waitAllTasks();
                auto ansCount = queryData->getCurrentAnswersCount();
                auto answers = queryData->getAnswers(0, std::min((size_t) 10, ansCount)).second;

                LOG(INFO) << "Found " << ansCount << " answers. Examples: " << [&answers]() {
                    std::stringstream s;
                    for (auto &ans: answers) {
                        s << ans << ' ';
                    }
                    return s.str();
                }();

                auto queryDuration = profilingTimer.stop();
                LOG(INFO) << queryDuration << " seconds spent to process molecule " << i * workers + worker << ": "
                          << *group[i].smiles;

                {
                    std::lock_guard<std::mutex> guard(mutex);
                    answerCounters[i * workers + worker] = ansCount;
                    times[i * workers + worker] = queryDuration;
                }
            }
        }

        inline void run() override {
            auto queries = loadQueriesFromFile(_inputFile, _fingerprintProvided, _searchData->getBaseLibrary());
            std::vector<float> times(queries.size());
            std::vector<size_t> answerCounters(queries.size());
            std::mutex mutex;
            size_t skipped = 0;

            std::vector<std::future<void>> tasks;
            auto queriesGrouped = splitQueries(std::move(queries), _workers);
            for (size_t worker = 0; worker < _workers; worker++) {
                auto &group = queriesGrouped[worker];
                tasks.emplace_back(std::async(std::launch::async, processGroup,
                                              std::cref(group), this->_searchData, std::ref(skipped),
                                              std::ref(mutex), worker, this->_workers, std::ref(answerCounters),
                                              std::ref(times)));
            }
            for (auto &task: tasks) {
                task.wait();
            }


            showStatistics(times, skipped, std::cout);
            showSummary(answerCounters, _summaryFile);
            ProfilingPool::showStatistics(std::cout);
        }


    };
} // qtr