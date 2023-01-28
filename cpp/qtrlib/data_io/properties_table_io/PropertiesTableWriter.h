#pragma once

#include "glog/logging.h"

#include "basic_io/BasicDataWriter.h"
#include "PropertiesFilter.h"

namespace qtr {

    class PropertiesTableWriter
            : public BasicDataWriter<std::pair<uint64_t, PropertiesFilter::property_list_t>, PropertiesTableWriter, std::ofstream> {
    private:
        uint64_t _writtenProperties;

    public:
        explicit PropertiesTableWriter(const std::filesystem::path &fileName) : BaseWriter(fileName),
                                                                                _writtenProperties(0) {
            _binaryWriter->write((char *) &_writtenProperties, sizeof _writtenProperties);
            LOG(INFO) << "Create properties table writer to " << fileName << " (" << _binaryWriter << ")";
        }

        ~PropertiesTableWriter() override {
            _binaryWriter->seekp(0, std::ios::beg);
            _binaryWriter->write((char *) &_writtenProperties, sizeof _writtenProperties);
            LOG(INFO) << "Delete properties table writer with " << _writtenProperties << " molecules (" << _binaryWriter
                      << ")";
        }

        PropertiesTableWriter &operator<<(const WriteValue &value) override {
            _writtenProperties++;
            auto &[id, propertiesList] = value;
            _binaryWriter->write((char *) &id, sizeof id);
            _binaryWriter->write((char *) &propertiesList, sizeof propertiesList);
            return *this;
        }

        using BaseWriter::operator<<;
    };

} // qtr
