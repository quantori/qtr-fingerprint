#include <utility>

#include "IndigoDriveFilter.h"
#include "string_table_io/StringTableReader.h"

using namespace std;

namespace qtr {
    string IndigoDriveFilter::getSmiles(CIDType id) {
        return _smilesTable.at(id);
    }

    void IndigoDriveFilter::initBallTreeLeaf(const filesystem::path &leafDirPath) {
        filesystem::path smilesTablePath = leafDirPath / ("smiles" + stringTableExtension);
        StringTableReader reader(smilesTablePath);
        _smilesTable = MapSmilesTable(reader.begin(), reader.end());
    }

    IndigoDriveFilter::IndigoDriveFilter(shared_ptr<const string> querySmiles) :
            IndigoFilter(std::move(querySmiles)) {}

    unique_ptr <ByIdAnswerFilter> IndigoDriveFilter::copy() {
        return make_unique<IndigoDriveFilter>(_querySmiles);
    }
} // qtr