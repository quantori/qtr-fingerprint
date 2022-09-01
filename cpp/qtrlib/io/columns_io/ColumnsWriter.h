#pragma once

#include "basic_io/BasicWriter.h"

namespace qtr {

    // TODO class is not tested after refactoring
    class ColumnsWriter : public BasicWriter<size_t, ColumnsWriter> {
    public:
        explicit ColumnsWriter(std::ostream *stream) : BaseWriter(stream) {}

        explicit ColumnsWriter(const std::filesystem::path &fileName) : ColumnsWriter(new std::ofstream(fileName)) {};

        void write(const size_t &value) override {
            if (_writtenSmiles != 0)
                *_stream << ' ';
            *_stream << value;
            _writtenSmiles++;
        }

        using BaseWriter::write;
    };

} // namespace qtr
