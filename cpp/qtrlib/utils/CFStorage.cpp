#include "CFStorage.h"

#include "Utils.h"

#include "base_cpp/scanner.h"
#include "molecule/cmf_loader.h"

using namespace qtr;
using namespace std;
using namespace indigo;

CFStorage::ValueType &CFStorage::Add(KeyType key, ValueType &&value) {
    size_t newSize = max((size_t)1, _map.size());
    while (key >= newSize)
        newSize <<= 1;
    if (_map.size() != newSize)
        _map.resize(newSize);
    _map[key] = make_unique<ValueType>(std::move(value));
    return *_map[key];
}

unique_ptr<Molecule> CFStorage::operator[](CFStorage::KeyType key) const {
    const ValueType &arr = *_map[key];
    BufferScanner scanner(arr);
    CmfLoader cmf_loader(scanner);
    auto molecule = make_unique<Molecule>();
    cmf_loader.loadMolecule(*molecule);
    return std::move(molecule);
}

CFStorage::CFStorage() : _map() {}
