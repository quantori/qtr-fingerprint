#include "FingerprintTableSearchEngine.h"

#include "QtrIndigoFingerprint.h"

#include "IndigoException.h"
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
    _serializedMolecules.clear();
    _fingerprintTable.clear();

    size_t moleculesNumber = 0;
    size_t failuresNumber = 0;
    IndigoSDFileIterator iterator = _indigoSessionPtr->iterateSDFile(path);

    for(IndigoMoleculeSPtr &molecule : iterator) {
        
        try {
            molecule->aromatize();

            byte *buf = nullptr; int size = 0;
            int result = indigoSerialize(molecule->id(), &buf, &size);
            
            _indigoSessionPtr->_checkResult(result);
            _serializedMolecules.emplace_back(buf, buf + size);

            QtrIndigoFingerprint fingerprint(*molecule, "sub");
            _fingerprintTable.emplace_back(fingerprint);
        }
        catch(const IndigoException &) {
            failuresNumber++;
        }

        moleculesNumber++;
        if (moleculesNumber % 1000 == 0)
            LOG(INFO) << "Processed " << moleculesNumber << " molecules...";
    }
    LOG(INFO) << "Processed " << moleculesNumber << " molecules (including " << failuresNumber << " failures)";
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
    LOG(INFO) << "Views for brute force: " << views.size();

    for(const IndigoFingerprintTableView *view : views) {
        for (IndigoFingerprintTableView::IndexType idx : *view) {
            
            qtr::IndigoFingerprint f = _fingerprintTable.at(idx);
            f &= fp;
            
            if (f.count() != bitsCount)
                continue;

            const std::vector<byte> &buf = _serializedMolecules.at(idx);
            const int moleculeId = indigoUnserialize(buf.data(), buf.size());

            _indigoSessionPtr->_checkResult(moleculeId);
            IndigoMolecule molecule(moleculeId, _indigoSessionPtr);

            IndigoSubstructureMatcher matcher = _indigoSessionPtr->substructureMatcher(molecule);
        
            if (!matcher.match(mol))
                continue;

            result.push_back(molecule);
        }
    }

    return result;
}

} // namespace qtr