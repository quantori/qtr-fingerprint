#pragma once

#include <utility>

#include "RunMode.h"
#include "iostream"
#include "search_data/QtrRamSearchData.h"
#include "search_data/QtrDriveSearchData.h"
#include "search_data/BingoNoSQLSearchData.h"

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
                SearchData::Query query(std::make_unique<std::string>(smiles), nullptr, _searchData->getBaseLibrary());
                ProfilingTimer profilingTimer("Query processing");
                try {
                    auto queryData = this->_searchData->search(query, PropertiesFilter::Bounds());
                    if (queryData == nullptr) {
                        LOG(ERROR) << "Can not parse given smiles";
                        continue;
                    }
                    queryData->waitAllTasks();
                    LOG(INFO) << "Found " << queryData->getCurrentAnswersCount() << " answers";
                    auto answers = queryData->getAnswers(0, 5).second;
                    LOG(INFO) << "Answer examples:";
                    if (dynamic_cast<const QtrDriveSearchData *>(this->_searchData.get()) != nullptr
                        || dynamic_cast<const QtrRamSearchData *>(this->_searchData.get()) != nullptr) {
                        for (auto &i: answers) {
                            LOG(INFO) << i;
                        }
                    } else if (dynamic_cast<const BingoNoSQLSearchData *>(this->_searchData.get()) != nullptr) {
                        for (auto &i: answers) {
                            LOG(INFO) << indigoSmiles(i);
                        }
                    }
                    std::cout << "Search time: " << profilingTimer.stop() << std::endl;
                } catch (std::exception &e) {
                    LOG(ERROR) << e.what() << "error occurred while processing " << smiles;
                    continue;
                }
            }
        }

    };
} // qtr