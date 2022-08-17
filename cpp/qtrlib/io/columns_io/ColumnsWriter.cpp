#include "ColumnsWriter.h"


namespace qtr {

    void ColumnsWriter::write(const std::vector<size_t> &columns) {
        for (size_t i = 0; i < columns.size(); i++) {
            *_outStream << columns[i];
            *_outStream << (i + 1 == columns.size() ? "" : " ");
        }
    }

    ColumnsWriter::~ColumnsWriter() {
        delete _outStream;
    }
} // namespace qtr
