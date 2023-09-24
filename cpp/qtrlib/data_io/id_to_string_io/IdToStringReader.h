#pragma once

#include "basic_io/BasicDataReader.h"
#include "io/BufferedReader.h"

namespace qtr {

    class IdToStringReader : public BasicDataReader<std::pair<uint64_t, std::string>, IdToStringReader, std::ifstream> {
    public:

        explicit IdToStringReader(const std::filesystem::path &fileName) : BaseReader(fileName) {}

        IdToStringReader &operator>>(ReadValue &value) override {
            auto& [id, str] = value;
            *_binaryReader >> id;
            char space = _binaryReader->get();
            assert(space == ' ');
            for (int symbol = _binaryReader->get(); symbol != EOF && symbol != '\n'; symbol = _binaryReader->get())  {
                str += (char)symbol;
            }
            return *this;
        }

        using BaseReader::operator>>;
    };

} // namespace qtr
