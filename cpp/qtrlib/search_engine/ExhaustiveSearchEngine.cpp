#include "ExhaustiveSearchEngine.h"
#include "QtrIndigoFingerprint.h"
#include "Utils.h"

#include "IndigoSDFileIterator.h"
#include "IndigoSubstructureMatcher.h"
#include "indigo.h"

#include <glog/logging.h>

using namespace indigo_cpp;

ExhaustiveSearchEngine::ExhaustiveSearchEngine(const IndigoSessionPtr &indigoSessionPtr)
    : _indigoSessionPtr(indigoSessionPtr)
{}

void ExhaustiveSearchEngine::build(const std::string &path)
{
    _moleculesInfo.clear();

    size_t moleculesNumber = 0;
    IndigoSDFileIterator iterator = _indigoSessionPtr->iterateSDFile(path);

    for(IndigoMoleculeSPtr &molecule : iterator) {
        molecule->aromatize();

        QtrIndigoFingerprint fingerprint(*molecule, "sub");

        MoleculeInfo moleculeInfo = {std::move(*molecule), std::move(fingerprint)};
        _moleculesInfo.push_back(std::move(moleculeInfo));

        moleculesNumber++;
        if (moleculesNumber % 1000 == 0)
            LOG(INFO) << "Processed " << moleculesNumber << " molecules...";
    }
}

std::vector<indigo_cpp::IndigoMolecule> ExhaustiveSearchEngine::findOverMolecules(const indigo_cpp::IndigoQueryMolecule &mol)
{
    std::vector<indigo_cpp::IndigoMolecule> result;

    indigoAromatize(mol.id());

    QtrIndigoFingerprint fingerprint(mol, "sub");
    int bitsCount = fingerprint.countBits();

    for(const MoleculeInfo &moleculeInfo : _moleculesInfo) {
        
        int commonBits = QtrIndigoFingerprint::commonBits(fingerprint, moleculeInfo.fingerprint);
        if (commonBits != bitsCount)
            continue;

        IndigoSubstructureMatcher matcher = _indigoSessionPtr->substructureMatcher(moleculeInfo.molecule);
        if (!matcher.match(mol))
            continue;

        result.push_back(moleculeInfo.molecule);
    }

    return result;
}

ExhaustiveSearchEngine::~ExhaustiveSearchEngine()
{
}
