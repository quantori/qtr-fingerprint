#pragma once

#include <utility>

#include "RunMode.h"
#include "iostream"
#include "RunDbUtils.h"

namespace qtr {
    class InteractiveMode : public RunMode {
        const qtr::BallTreeSearchEngine &_ballTree;
        std::shared_ptr<const SmilesTable> _smilesTable;
        qtr::TimeTicker &_timeTicker;
        uint64_t _ansCount;
        uint64_t _threadsCount;
    public:
        inline InteractiveMode(const qtr::BallTreeSearchEngine &ballTree, shared_ptr<const SmilesTable> smilesTable,
                               qtr::TimeTicker &timeTicker, uint64_t ansCount, uint64_t threadsCount) :
                _ballTree(ballTree),
                _smilesTable(std::move(smilesTable)),
                _timeTicker(timeTicker),
                _ansCount(ansCount),
                _threadsCount(threadsCount) {}


        inline void run() override {
            while (true) {
                std::cout << "Enter smiles: ";
                std::string smiles;
                std::cin >> smiles;
                if (smiles.empty())
                    break;
                _timeTicker.tick();
                try {
                    auto [error, queryData] = doSearch(smiles, _ballTree, _smilesTable, _ansCount, _threadsCount);
                    if (error) {
                        LOG(ERROR) << "Can not parse given smiles";
                        continue;
                    }
                    queryData->waitAllTasks();
                    LOG(INFO) << "Found " << queryData->getCurrentAnswersCount() << " answers";
                    auto answers = queryData->getAnswers(0, 5).second;
                    for (auto& i : answers) {
                        LOG(INFO) << (*_smilesTable)[i];
                    }
                    std::cout << "Search time: " << _timeTicker.tick("Search time") << std::endl;
                } catch (std::exception &e) {
                    LOG(ERROR) << e.what() << " while processing " << smiles;
                    continue;
                }
            }
        }

    };
} // qtr