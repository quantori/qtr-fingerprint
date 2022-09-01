#pragma once

#include "basic_io/BasicDataReader.h"

namespace qtr {

    // TODO class is not tested after refactoring
    class ColumnsReader : public BasicDataReader<size_t, ColumnsReader> {
    public:
        explicit ColumnsReader(std::istream *stream) : BaseReader(stream) {};

        explicit ColumnsReader(const std::filesystem::path &fileName) : ColumnsReader(new std::ifstream(fileName)) {}

        size_t readOne() override {
            size_t result;
            *_stream >> result;
            return result;
        }

        bool isEof() const override {
            return _stream->eof();
        }
    };

} // namespace qtr
