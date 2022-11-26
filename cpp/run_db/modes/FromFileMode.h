#pragma once

#include "RunMode.h"
#include "RunDbUtils.h"

namespace qtr {
    class FromFileMode : public RunMode {
        const qtr::BallTreeSearchEngine &ballTree;
        SmilesTable &smilesTable;
        qtr::TimeTicker &timeTicker;
        const std::filesystem::path &inputFile;
        uint64_t ansCount;
        uint64_t startSearchDepth;
    public:
        inline FromFileMode(const qtr::BallTreeSearchEngine &ballTree, SmilesTable &smilesTable,
                            qtr::TimeTicker &timeTicker, const std::filesystem::path &inputFile, uint64_t ansCount,
                            uint64_t startSearchDepth) :
                ballTree(ballTree),
                smilesTable(smilesTable),
                timeTicker(timeTicker),
                inputFile(inputFile),
                ansCount(ansCount),
                startSearchDepth(startSearchDepth) {}


        inline void run() override {
            std::ifstream input(inputFile);
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
                timeTicker.tick();
                const auto result = doSearch(queries[i], ballTree, smilesTable, ansCount, startSearchDepth);
                if (result.first) {
                    skipped++;
                    continue;
                }
                times.emplace_back(timeTicker.tick("search molecule " + std::to_string(i) + ": " + queries[i]));
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
        }

    };
} // qtr