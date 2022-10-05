#include "SmilesRandomAccessTable.h"
#include "indexes_io/IndexesWriter.h"
#include "indexes_io/IndexesReader.h"
#include "smiles_table_io/SmilesTableReader.h"

namespace qtr {

    namespace qtr {
        void buildSmilesFile(const std::filesystem::path &smilesTablePath, const std::filesystem::path &tableDir,
                             std::vector<uint64_t> &seeks) {
            uint64_t writtenNumber = 0;
            uint64_t writtenBytes = 0;
            std::ofstream smilesOut(SmilesRandomAccessTable::getSmilesStoragePath(tableDir));
            for (const auto &[id, smiles]: SmilesTableReader(smilesTablePath)) {
                assert(writtenNumber == id);
                writtenBytes += smiles.size();
                writtenNumber++;
                seeks.emplace_back(writtenBytes);
                smilesOut.write(smiles.c_str(), smiles.size());
            }
        }
    }

    SmilesRandomAccessTable::SmilesRandomAccessTable(const std::filesystem::path &tableDir) {
        _smilesFile = std::ifstream(SmilesRandomAccessTable::getSmilesStoragePath(tableDir));
        IndexesReader seeksIn(SmilesRandomAccessTable::getSeeksPath(tableDir));
        seeksIn >> _seeks;
    }

    std::string SmilesRandomAccessTable::operator[](size_t i) {
        std::string result;
        uint64_t begin = i == 0 ? 0 : _seeks[i - 1];
        uint64_t end = _seeks[i];
        result.resize(end - begin);
        _smilesFile.seekg(begin);
        _smilesFile.read(const_cast<char *>(result.c_str()), end - begin);
        return result;
    }

    SmilesRandomAccessTable::SmilesRandomAccessTable(const std::filesystem::path &smilesTablePath,
                                                     const std::filesystem::path &tableDir) {
        qtr::buildSmilesFile(smilesTablePath, tableDir, _seeks);
        IndexesWriter seeksOut(SmilesRandomAccessTable::getSeeksPath(tableDir));
        seeksOut << _seeks;
        _smilesFile = std::ifstream(SmilesRandomAccessTable::getSmilesStoragePath(tableDir));
    }
}