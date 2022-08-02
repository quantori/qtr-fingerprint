#include "ColumnsWriter.h"


namespace qtr {

    void ColumnsWriter::write(const std::vector<column_t> &columns) {
        for (size_t i = 0; i < columns.size(); i++) {
            *_outStream << columns[i];
            *_outStream << (i + 1 == columns.size() ? "" : " ");
        }
    }

    ColumnsWriter::~ColumnsWriter() {
        delete _outStream;
    }

    ColumnsWriter::ColumnsWriter(std::filesystem::path fileName) {
        fileName.replace_extension(columnsExtension);
        _outStream = new std::ofstream(fileName);
    }

} // namespace qtr
