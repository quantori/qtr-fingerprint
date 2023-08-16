#pragma once

#include <string>

#include "glog/logging.h"

#include "StringTableIOConsts.h"
#include "basic_io/BasicDataWriter.h"

namespace qtr {

    class StringTableWriter
            : public BasicDataWriter<string_table_value_t, StringTableWriter, std::ofstream> {
    private:
        uint64_t _writtenStrings;

    public:
        explicit StringTableWriter(const std::filesystem::path &fileName) : BaseWriter(fileName), _writtenStrings(0) {
            _binaryWriter->write((char *) &_writtenStrings, sizeof _writtenStrings); // reserve space for table size
            LOG(INFO) << "Create string table writer to " << fileName << " (" << _binaryWriter << ")";
        }

        ~StringTableWriter() override {
            _binaryWriter->seekp(0, std::ios::beg);
            _binaryWriter->write((char *) &_writtenStrings, sizeof _writtenStrings);
            LOG(INFO) << "Delete string table writer with " << _writtenStrings << " molecules (" << _binaryWriter
                      << ")";
        }

        StringTableWriter &operator<<(const WriteValue &value) override {
            _writtenStrings++;
            auto &[id, str] = value;
            _binaryWriter->write((char *) &id, sizeof id);
            _binaryWriter->write(str.c_str(), str.size());
            _binaryWriter->write("\n", 1);
            return *this;
        }

        using BaseWriter::operator<<;
    };

} // qtr