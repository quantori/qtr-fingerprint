#include "DecisionTreeSearchEngine.h"

#include "QtrIndigoFingerprint.h"

#include "IndigoSDFileIterator.h"

#include <glog/logging.h>

using namespace indigo_cpp;

namespace qtr{

DecisionTreeSearchEngine::DecisionTreeSearchEngine(const indigo_cpp::IndigoSessionPtr &indigoSessionPtr)
    : _indigoSessionPtr(indigoSessionPtr)
{}

DecisionTreeSearchEngine::~DecisionTreeSearchEngine()
{}

void DecisionTreeSearchEngine::build(const std::string &path)
{
    _molecules.clear();
    _fingerprintTable.clear();

    size_t moleculesNumber = 0;
    IndigoSDFileIterator iterator = _indigoSessionPtr->iterateSDFile(path);

    for(IndigoMoleculeSPtr &molecule : iterator) {
        
        molecule->aromatize();
        _molecules.push_back(std::move(*molecule));

        QtrIndigoFingerprint fingerprint(*molecule, "sub");
        _fingerprintTable.push_back(qtr::IndigoFingerprint());
        _fingerprintTable.back().setBytes(fingerprint.data());

        moleculesNumber++;
        if (moleculesNumber % 1000 == 0)
            LOG(INFO) << "Processed " << moleculesNumber << " molecules...";
    }
}


} // namespace qtr