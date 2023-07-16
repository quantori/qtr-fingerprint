#pragma once

#include <utility>

#include "RunMode.h"
#include "iostream"
#include "search_data/QtrRamSearchData.h"
#include "search_data/QtrDriveSearchData.h"

namespace qtr {
    class InteractiveMode : public RunMode {
    public:
        inline explicit InteractiveMode(std::shared_ptr<SearchData> searchData) :
                RunMode(std::move(searchData)) {}


        inline void run() override {
            while (true) {
                std::cout << "Enter smiles: ";
                std::string smiles;
                std::cin >> smiles;
                if (smiles.empty())
                    break;
                this->_searchData->timeTicker.tick();
                try {
                    auto queryData = this->_searchData->search(smiles, PropertiesFilter::Bounds());
                    if (queryData == nullptr) {
                        LOG(ERROR) << "Can not parse given smiles";
                        continue;
                    }
                    queryData->waitAllTasks();
                    LOG(INFO) << "Found " << queryData->getCurrentAnswersCount() << " answers";
                    auto answers = queryData->getAnswers(0, 5).second;
                    if (dynamic_cast<const QtrRamSearchData *>(this->_searchData.get()) != nullptr) {
                        LOG(INFO) << "Answer examples:";
                        const auto *ramSearchData = dynamic_cast<const QtrRamSearchData *>(this->_searchData.get());
                        for (auto &i: answers) {
                            LOG(INFO) << (*ramSearchData->smilesTable)[i];
                        }
                    } else if (dynamic_cast<const QtrDriveSearchData *>(this->_searchData.get()) != nullptr) {
                        LOG(INFO) << "Answer examples:";
                        for (auto &i: answers) {
                            LOG(INFO) << i;
                        }
                    }
                    std::cout << "Search time: " << this->_searchData->timeTicker.tick("Search time") << std::endl;
                } catch (std::exception &e) {
                    LOG(ERROR) << e.what() << " while processing " << smiles;
                    continue;
                }
            }
        }

    };
} // qtr