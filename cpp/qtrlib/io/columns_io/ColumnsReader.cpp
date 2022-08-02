#include "ColumnsReader.h"

namespace qtr {

    std::vector<column_t> qtr::ColumnsReader::read() {
        std::vector<column_t> columns;
        column_t col;
        while (*_inStream >> col) {
            columns.emplace_back(col);
        }
        return columns;
    }

    ColumnsReader::~ColumnsReader() {
        delete _inStream;
    }

} // namespace qtr
