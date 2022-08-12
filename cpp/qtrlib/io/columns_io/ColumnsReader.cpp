#include "ColumnsReader.h"

namespace qtr {

    std::vector<size_t> qtr::ColumnsReader::readAll() {
        std::vector<size_t> columns;
        size_t col;
        while (*_inStream >> col) {
            columns.emplace_back(col);
        }
        return columns;
    }

    ColumnsReader::~ColumnsReader() {
        delete _inStream;
    }

} // namespace qtr
