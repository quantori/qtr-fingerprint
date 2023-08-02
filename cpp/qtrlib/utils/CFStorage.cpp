#include "CFStorage.h"

#include "Utils.h"

#include "base_cpp/scanner.h"
#include "molecule/cmf_loader.h"

using namespace qtr;
using namespace std;
using namespace indigo;

CFStorage::ValueType &CFStorage::Add(KeyType key, ValueType &&value) {
    if (key != _map.size()) {
        logErrorAndExit("CFStorage: keys should be added in order 0, 1, 2, 3, ...");
        // TODO: avoid this restriction
    }
    _map.emplace_back(std::move(value));
    return _map.back();
}

unique_ptr<Molecule> CFStorage::operator[](CFStorage::KeyType key) const {
    const ValueType &arr = _map[key];
    BufferScanner scanner(arr);
    CmfLoader cmf_loader(scanner);
    unique_ptr<Molecule> molecule;
    cmf_loader.loadMolecule(*molecule);
    return std::move(molecule);
}

CFStorage::CFStorage() : _map() {}
