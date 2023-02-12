#pragma once

#include <string>

#include "StringTableIOConsts.h"
#include "basic_io/BasicDataReader.h"
#include "io/BufferedReader.h"

namespace qtr {

    class StringTableReader
            : public BasicDataReader<string_table_value_t, StringTableReader, BufferedReader<>> {
    private:
        uint64_t _inStream;

    public:
        explicit StringTableReader(const std::filesystem::path &fileName)
                : BaseReader(fileName), _inStream(0) {
            _binaryReader->read((char *) &_inStream, sizeof _inStream);
            LOG(INFO) << "Create string table reader with " << _inStream << " string (" << _binaryReader << ")";
        }

        ~StringTableReader() override {
            LOG(INFO) << "Delete string table reader (" << _binaryReader << ")";
        }

        StringTableReader &operator>>(ReadValue &value) override {
            auto &[id, str] = value;
            _binaryReader->read((char *) &id, sizeof id);
            char symbol;
            while ((symbol = (char) _binaryReader->get()) != '\n') {
                str += symbol;
            }
            return *this;
        }

        using BaseReader::operator>>;
    };


} // qtr