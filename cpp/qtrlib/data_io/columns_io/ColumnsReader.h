#pragma once

#include "basic_io/BasicDataReader.h"
#include "io/BufferedReader.h"

namespace qtr {

    // TODO class is not tested after refactoring
    class ColumnsReader : public BasicDataReader<size_t, ColumnsReader, std::ifstream> {
    public:

        explicit ColumnsReader(const std::filesystem::path &fileName) : BaseReader(fileName) {}

        size_t readOne() override {
            size_t result;
            *_binaryReader >> result;
            return result;
        }

        bool eof() const override {
            return _binaryReader->eof();
        }
    };

} // namespace qtr
