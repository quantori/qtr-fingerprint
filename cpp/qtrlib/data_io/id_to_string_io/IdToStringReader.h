#pragma once

#include "basic_io/BasicDataReader.h"
#include "io/BufferedReader.h"

namespace qtr {

    class IdToStringReader : public BasicDataReader<std::pair<uint64_t, std::string>, IdToStringReader, std::ifstream> {
    public:

        explicit IdToStringReader(const std::filesystem::path &fileName) : BaseReader(fileName) {}

        IdToStringReader &operator>>(ReadValue &value) override {
            *_binaryReader >> value.first >> value.second;
            return *this;
        }

        using BaseReader::operator>>;
    };

} // namespace qtr
