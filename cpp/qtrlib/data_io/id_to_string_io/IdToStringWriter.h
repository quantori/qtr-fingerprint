#pragma once

#include "basic_io/BasicDataWriter.h"

namespace qtr {

    class IdToStringWriter : public BasicDataWriter<std::pair<uint64_t, std::string>, IdToStringWriter, std::ofstream> {
    public:
        explicit IdToStringWriter(const std::filesystem::path &fileName) : BaseWriter(fileName), _writtenRows(0) {}

        IdToStringWriter &operator<<(const WriteValue &value) override {
            if (_writtenRows != 0)
                *_binaryWriter << '\n';
            *_binaryWriter << value.first << ' ' << value.second;
            _writtenRows++;
            return *this;
        }

        using BaseWriter::operator<<;

    private:
        size_t _writtenRows;
    };

} // namespace qtr