#pragma once

#include <filesystem>

#include "IndigoFilter.h"
#include "MapSmilesTable.h"
#include "QtrBallTreeLeafInitMixin.h"

namespace qtr {

    class IndigoDriveFilter : public IndigoFilter, public QtrBallTreeLeafInitMixin {
    public:
        void initBallTreeLeaf(const std::filesystem::path &leafDirPath) override;

        std::unique_ptr<AnswerFilter> copy() override;

        explicit IndigoDriveFilter(std::shared_ptr<const std::string> querySmiles);

    private:
        std::string getSmiles(CIDType id) override;

        MapSmilesTable _smilesTable;
    };

} // qtr

