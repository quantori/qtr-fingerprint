#pragma once

#include "basic_io/BasicDataWriter.h"

namespace qtr {

    class ColumnsWriter : public BasicDataWriter<size_t, ColumnsWriter, std::ofstream> {
    public:

        explicit ColumnsWriter(const std::filesystem::path &fileName) : BaseWriter(fileName), _writtenColumns(0) {};

        ColumnsWriter &operator<<(const WriteValue &value) override {
            if (_writtenColumns != 0)
                *_binaryWriter << ' ';
            *_binaryWriter << value;
            _writtenColumns++;
            return *this;
        }

        using BaseWriter::operator<<;

        size_t _writtenColumns;
    private:
    };

} // namespace qtr
