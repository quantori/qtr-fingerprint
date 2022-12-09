#pragma once

#include "basic_io/BasicDataWriter.h"

namespace qtr {

    class IndexesWriter : public BasicDataWriter<uint64_t, IndexesWriter, std::ofstream> {
    public:

        explicit IndexesWriter(const std::filesystem::path &fileName) : BaseWriter(fileName), _writtenColumns(0) {};

        IndexesWriter &operator<<(const WriteValue &value) override {
            if (_writtenColumns != 0)
                *_binaryWriter << ' ';
            *_binaryWriter << value;
            _writtenColumns++;
            return *this;
        }

        using BaseWriter::operator<<;

    private:
        size_t _writtenColumns;
    };

} // namespace qtr
