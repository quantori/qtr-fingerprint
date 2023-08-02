#include <utility>

#include "IndigoSmilesDriveFilter.h"
#include "string_table_io/StringTableReader.h"

using namespace std;

namespace qtr {
    string IndigoSmilesDriveFilter::getSmiles(CIDType id) {
        return _smilesTable.at(id);
    }

    void IndigoSmilesDriveFilter::initBallTreeLeaf(const filesystem::path &leafDirPath) {
        filesystem::path smilesTablePath = leafDirPath / ("smiles" + stringTableExtension);
        StringTableReader reader(smilesTablePath);
        _smilesTable = MapSmilesTable(reader.begin(), reader.end());
    }

    IndigoSmilesDriveFilter::IndigoSmilesDriveFilter(shared_ptr<const string> querySmiles) :
            IndigoSmilesFilter(std::move(querySmiles)) {}

    unique_ptr <ByIdAnswerFilter> IndigoSmilesDriveFilter::copy() {
        return make_unique<IndigoSmilesDriveFilter>(_querySmiles);
    }
} // qtr