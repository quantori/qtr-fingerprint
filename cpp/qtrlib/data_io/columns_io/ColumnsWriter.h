#pragma once

#include "basic_io/BasicDataWriter.h"

namespace qtr {

    class ColumnsWriter : public BasicDataWriter<size_t, ColumnsWriter, std::ofstream> {
    public:

        explicit ColumnsWriter(const std::filesystem::path &fileName) : BaseWriter(fileName), _writtenColumns(0) {};

        void write(const size_t &value) override {
            if (_writtenColumns != 0)
                *_binaryWriter << ' ';
            *_binaryWriter << value;
            _writtenColumns++;
        }

        using BaseWriter::write;

        size_t _writtenColumns;
    private:
    };

} // namespace qtr
