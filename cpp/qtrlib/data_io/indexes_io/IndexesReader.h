#pragma once

#include "basic_io/BasicDataReader.h"
#include "io/BufferedReader.h"

namespace qtr {

    class IndexesReader : public BasicDataReader<uint64_t, IndexesReader, std::ifstream> {
    public:

        explicit IndexesReader(const std::filesystem::path &fileName) : BaseReader(fileName) {}

        IndexesReader &operator>>(ReadValue &value) override {
            *_binaryReader >> value;
            return *this;
        }

        using BaseReader::operator>>;
    };

} // namespace qtr
