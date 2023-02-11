#pragma once

#include "glog/logging.h"

#include "io/BufferedReader.h"
#include "basic_io/BasicDataReader.h"
#include "PropertiesFilter.h"

namespace qtr {

    class PropertiesTableReader :
            public BasicDataReader<std::pair<uint64_t, PropertiesFilter::Properties>, PropertiesTableReader, BufferedReader<>> {
    private:
        uint64_t _propertiesInStream;

    public:
        explicit PropertiesTableReader(const std::filesystem::path &fileName) : BaseReader(fileName),
                                                                                _propertiesInStream(0) {
            _binaryReader->read((char *) &_propertiesInStream, sizeof _propertiesInStream);
            LOG(INFO) << "Create properties table reader with " << _propertiesInStream << " SMILES (" << _binaryReader
                      << ")";
        }

        ~PropertiesTableReader() override {
            LOG(INFO) << "Delete properties table reader (" << _binaryReader << ")";
        }

        PropertiesTableReader &operator>>(ReadValue &value) override {
            auto &[id, propertiesList] = value;
            _binaryReader->read((char *) &id, sizeof id);
            _binaryReader->read((char *) &propertiesList, sizeof propertiesList);
            return *this;
        }

        using BaseReader::operator>>;
    };
}
