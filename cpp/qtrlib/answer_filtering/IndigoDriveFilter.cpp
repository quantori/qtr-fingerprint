#include <utility>
#include <memory>

#include "IndigoDriveFilter.h"
#include "string_table_io/StringTableReader.h"
#include "IndigoMolecule.h"
#include "indigo.h"


using namespace std;

namespace qtr {
    void IndigoDriveFilter::initBallTreeLeaf(const filesystem::path &leafDirPath) {
        filesystem::path smilesTablePath = leafDirPath / ("smiles" + stringTableExtension);
        StringTableReader reader(smilesTablePath);
        _smilesTable = MapSmilesTable(reader.begin(), reader.end());
    }

    IndigoDriveFilter::IndigoDriveFilter(shared_ptr<const string> querySmiles) :
            IndigoFilter(std::move(querySmiles)) {}

    unique_ptr <AnswerFilter> IndigoDriveFilter::copy() {
        return make_unique<IndigoDriveFilter>(_querySmiles);
    }

    std::shared_ptr<indigo_cpp::IndigoMolecule> IndigoDriveFilter::getMolecule(CIDType id) {
        auto mol = std::make_shared<indigo_cpp::IndigoMolecule>(_indigoSessionPtr->loadMolecule(_smilesTable.at(id)));
        mol->aromatize();
        return mol;
    }
} // qtr