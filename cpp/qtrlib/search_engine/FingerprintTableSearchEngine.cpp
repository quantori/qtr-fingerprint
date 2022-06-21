#include "FingerprintTableSearchEngine.h"

#include "QtrIndigoFingerprint.h"

#include "IndigoSDFileIterator.h"
#include "IndigoSubstructureMatcher.h"

#include "indigo.h"

#include <glog/logging.h>

using namespace indigo_cpp;

namespace qtr {

FingerprintTableSearchEngine::FingerprintTableSearchEngine(const indigo_cpp::IndigoSessionPtr &indigoSessionPtr)
    : _indigoSessionPtr(indigoSessionPtr)
{}

FingerprintTableSearchEngine::~FingerprintTableSearchEngine()
{}

void FingerprintTableSearchEngine::build(const std::string &path)
{
    _molecules.clear();
    _fingerprintTable.clear();

    size_t moleculesNumber = 0;
    IndigoSDFileIterator iterator = _indigoSessionPtr->iterateSDFile(path);

    for(IndigoMoleculeSPtr &molecule : iterator) {
        
        molecule->aromatize();
        _molecules.push_back(std::move(*molecule));

        QtrIndigoFingerprint fingerprint(_molecules.back(), "sub");
        _fingerprintTable.push_back(qtr::IndigoFingerprint());
        _fingerprintTable.back().setBytes(fingerprint.data());

        moleculesNumber++;
        if (moleculesNumber % 1000 == 0)
            LOG(INFO) << "Processed " << moleculesNumber << " molecules...";
    }
}

std::vector<indigo_cpp::IndigoMolecule> FingerprintTableSearchEngine::findOverMolecules(const indigo_cpp::IndigoQueryMolecule &mol)
{
    std::vector<IndigoMolecule> result;

    indigoAromatize(mol.id());

    QtrIndigoFingerprint fingerprint(mol, "sub");
    int bitsCount = fingerprint.countBits();

    qtr::IndigoFingerprint fp;
    fp.setBytes(fingerprint.data());

    std::vector<const IndigoFingerprintTableView *> views = findTableViews(fp);

    for(const IndigoFingerprintTableView *view : views) {
        for (IndigoFingerprintTableView::IndexType idx : *view) {
            
            qtr::IndigoFingerprint f = _fingerprintTable.at(idx);
            f &= fp;
            
            if (f.count() != bitsCount)
                continue;
            
            const IndigoMolecule &molecule = _molecules.at(idx);
            IndigoSubstructureMatcher matcher = _indigoSessionPtr->substructureMatcher(molecule);
        
            if (!matcher.match(mol))
                continue;

            result.push_back(molecule);
        }
    }

    return result;
}

} // namespace qtr