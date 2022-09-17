#pragma once

#include "basic_io/BasicDataReader.h"
#include "FingerprintTableIOConsts.h"


namespace qtr {

    // TODO: test this class
    class FingerprintTableReader
            : public BasicDataReader<fingerprint_table_value_t, FingerprintTableReader, std::ifstream> {

    private:
        uint64_t _fingerprintsInStream;

    public:
        explicit FingerprintTableReader(const std::filesystem::path &fileName) : BaseReader(fileName),
                                                                                 _fingerprintsInStream(0) {
            _binaryReader->read((char *) &_fingerprintsInStream, sizeof _fingerprintsInStream);
            LOG(INFO) << "Create fingerprint table reader with " << _fingerprintsInStream << " SMILES ("
                      << _binaryReader << ")";
        }

        ~FingerprintTableReader() override {
            LOG(INFO) << "Delete fingerprint table reader (" << _binaryReader << ")";
        }

        FingerprintTableReader &operator>>(ReadValue &readValue) override {
            assert(_fingerprintsInStream > 0);
            auto &[id, fingerprint] = readValue;
            _binaryReader->read((char *) &id, sizeof id);
            fingerprint.load(*_binaryReader);
            _fingerprintsInStream--;
            return *this;
        }

        using BaseReader::operator>>;
    };

} // qtr