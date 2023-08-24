#pragma once

#include "RunMode.h"
#include "web/WebMode.h"
#include "FromFileMode.h"
#include "SearchData.h"
#include "RunArgs.h"
#include "InteractiveMode.h"

namespace qtr {

    class RunModeFactory {
    public:
        static std::unique_ptr<RunMode> create(const RunArgs &args, const std::shared_ptr<SearchData> &searchData) {
            std::unique_ptr<RunMode> result = nullptr;
            if (args.mode() == RunMode::Type::Interactive)
                result = make_unique<InteractiveMode>(searchData);
            else if (args.mode() == RunMode::Type::FromFile)
                result = make_unique<FromFileMode>(searchData, args.queriesFile());
            else if (args.mode() == RunMode::Type::Web)
                result = make_unique<WebMode>(searchData);
            return result;
        }
    };

} // qtr