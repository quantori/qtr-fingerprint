#pragma once

#include <utility>

#include "RunMode.h"
#include "iostream"
#include "RunDbUtils.h"

namespace qtr {
    class InteractiveMode : public RunMode {
        const qtr::BallTreeSearchEngine &ballTree;
        std::shared_ptr<const SmilesTable> smilesTable;
        qtr::TimeTicker &timeTicker;
        uint64_t ansCount;
        uint64_t startSearchDepth;
    public:
        inline InteractiveMode(const qtr::BallTreeSearchEngine &ballTree, shared_ptr<const SmilesTable> smilesTable,
                               qtr::TimeTicker &timeTicker, uint64_t ansCount, uint64_t startSearchDepth) :
                ballTree(ballTree),
                smilesTable(std::move(smilesTable)),
                timeTicker(timeTicker),
                ansCount(ansCount),
                startSearchDepth(startSearchDepth) {}


        inline void run() override {
            while (true) {
                std::cout << "Enter smiles: ";
                std::string smiles;
                std::cin >> smiles;
                if (smiles.empty())
                    break;
                timeTicker.tick();
                try {
                    auto [error, queryData] = doSearch(smiles, ballTree, smilesTable, startSearchDepth);
                    if (error) {
                        LOG(ERROR) << "Can not parse given smiles";
                        continue;
                    }
                    queryData->waitAllTasks();
                    LOG(INFO) << "Found " << queryData->getCurrentAnswersCount() << " answers";
                    auto answers = queryData->getAnswers(0, 5).second;
                    for (auto& i : answers) {
                        LOG(INFO) << (*smilesTable)[i];
                    }
                    std::cout << "Search time: " << timeTicker.tick("Search time") << std::endl;
                } catch (std::exception &e) {
                    LOG(ERROR) << e.what() << " while processing " << smiles;
                    continue;
                }
            }
        }

    };
} // qtr