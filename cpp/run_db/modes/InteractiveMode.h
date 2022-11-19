#pragma once

#include "RunMode.h"
#include "iostream"
#include "RunDbUtils.h"

namespace qtr {
    class InteractiveMode : public RunMode {
        const qtr::BallTreeSearchEngine &ballTree;
        SmilesTable &smilesTable;
        qtr::TimeTicker &timeTicker;
        uint64_t ansCount;
        uint64_t startSearchDepth;
    public:
        inline InteractiveMode(const qtr::BallTreeSearchEngine &ballTree, SmilesTable &smilesTable,
                               qtr::TimeTicker &timeTicker, uint64_t ansCount, uint64_t startSearchDepth) :
                ballTree(ballTree),
                smilesTable(smilesTable),
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
                const auto result = doSearch(smiles, ballTree, smilesTable, ansCount, startSearchDepth);
                if (result.first)
                    continue;
                std::cout << "Search time: " << timeTicker.tick("Search time") << std::endl;
            }
        }

    };
} // qtr