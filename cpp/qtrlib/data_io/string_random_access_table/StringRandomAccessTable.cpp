#include "StringRandomAccessTable.h"
#include "indexes_io/IndexesWriter.h"
#include "indexes_io/IndexesReader.h"
#include "string_table_io/StringTableReader.h"

namespace qtr {

    namespace qtr {
        void buildStringsFile(const std::filesystem::path &tablePath, const std::filesystem::path &tableDir,
                              std::vector<uint64_t> &seeks) {
            uint64_t writtenNumber = 0;
            uint64_t writtenBytes = 0;
            std::ofstream out(StringRandomAccessTable::getStoragePath(tableDir));
            for (const auto &[id, str]: StringTableReader(tablePath)) {
                assert(writtenNumber == id);
                writtenBytes += str.size();
                writtenNumber++;
                seeks.emplace_back(writtenBytes);
                out.write(str.c_str(), str.size());
            }
        }
    }

    StringRandomAccessTable::StringRandomAccessTable(const std::filesystem::path &tableDir) {
        _tableFile = std::ifstream(StringRandomAccessTable::getStoragePath(tableDir));
        IndexesReader seeksIn(StringRandomAccessTable::getSeeksPath(tableDir));
        seeksIn >> _seeks;
    }

    std::string StringRandomAccessTable::operator[](size_t i) {
        std::string result;
        uint64_t begin = i == 0 ? 0 : _seeks[i - 1];
        uint64_t end = _seeks[i];
        result.resize(end - begin);
        _tableFile.seekg(begin);
        _tableFile.read(const_cast<char *>(result.c_str()), end - begin);
        return result;
    }

    StringRandomAccessTable::StringRandomAccessTable(const std::filesystem::path &tablePath,
                                                     const std::filesystem::path &tableDir) {
        qtr::buildStringsFile(tablePath, tableDir, _seeks);
        IndexesWriter seeksOut(StringRandomAccessTable::getSeeksPath(tableDir));
        seeksOut << _seeks;
        _tableFile = std::ifstream(StringRandomAccessTable::getStoragePath(tableDir));
    }
}