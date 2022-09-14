#pragma once

#include "basic_io/BasicDataReader.h"
#include "io/BufferedReader.h"

namespace qtr {

    class ColumnsReader : public BasicDataReader<size_t, ColumnsReader, std::ifstream> {
    public:

        explicit ColumnsReader(const std::filesystem::path &fileName) : BaseReader(fileName) {}

        ColumnsReader &operator>>(size_t &value) override {
            *_binaryReader >> value;
            return *this;
        }

        using BaseReader::operator>>;
    };

} // namespace qtr
