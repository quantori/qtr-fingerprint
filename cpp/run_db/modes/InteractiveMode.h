#pragma once

#include <utility>

#include "RunMode.h"
#include "iostream"
#include "RunDbUtils.h"

namespace qtr {
    class InteractiveMode : public RunMode {
    private:
        std::shared_ptr<const SearchData> _searchData;

    public:
        inline explicit InteractiveMode(std::shared_ptr<const SearchData> searchData)
                : _searchData(std::move(searchData)) {}


        inline void run() override {
            while (true) {
                std::cout << "Enter smiles: ";
                std::string smiles;
                std::cin >> smiles;
                if (smiles.empty())
                    break;
                _searchData->timeTicker.tick();
                try {
                    auto [error, queryData] = runSearch(*_searchData, smiles, PropertiesFilter::Bounds());
                    if (error) {
                        LOG(ERROR) << "Can not parse given smiles";
                        continue;
                    }
                    queryData->waitAllTasks();
                    LOG(INFO) << "Found " << queryData->getCurrentAnswersCount() << " answers";
                    auto answers = queryData->getAnswers(0, 5).second;
                    if (_searchData->getClass() == SearchData::DerivedClasses::RamSearchData) {
                        const auto *ramSearchData = dynamic_cast<const RamSearchData *>(_searchData.get());
                        for (auto &i: answers) {
                            LOG(INFO) << (*ramSearchData->smilesTable)[i];
                        }
                    }
                    std::cout << "Search time: " << _searchData->timeTicker.tick("Search time") << std::endl;
                } catch (std::exception &e) {
                    LOG(ERROR) << e.what() << " while processing " << smiles;
                    continue;
                }
            }
        }

    };
} // qtr