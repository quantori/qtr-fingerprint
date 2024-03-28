#pragma once

#include "basic_io/BasicDataReader.h"
#include "FingerprintTableIOConsts.h"
#include "io/BufferedReader.h"


namespace qtr {

    class FingerprintTableReader
            : public BasicDataReader<fingerprint_table_value_t, FingerprintTableReader, BufferedReader<>> {

    private:
        uint64_t _fingerprintsInStream;
        uint64_t _fingerprintSize;

    public:
        explicit FingerprintTableReader(const std::filesystem::path &fileName) : BaseReader(fileName),
                                                                                 _fingerprintsInStream(0),
                                                                                 _fingerprintSize(0) {
            _binaryReader->read((char *) &_fingerprintsInStream, sizeof _fingerprintsInStream);
            _binaryReader->read((char *) &_fingerprintSize, sizeof _fingerprintSize);
        }

        FingerprintTableReader &operator>>(ReadValue &readValue) override {
            assert(_fingerprintsInStream > 0);
            auto &[id, fingerprint] = readValue;
            _binaryReader->read((char *) &id, sizeof id);
            fingerprint.load(*_binaryReader, _fingerprintSize);
            _fingerprintsInStream--;
            return *this;
        }

        [[nodiscard]] uint64_t fingerprintSize() const {
            return _fingerprintSize;
        }

        using BaseReader::operator>>;
    };

} // qtr