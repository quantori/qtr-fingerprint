#pragma once

#include "basic_io/BasicDataReader.h"
#include "Fingerprint.h"

namespace qtr {

    // TODO: test this class
    class FingerprintTableReader
            : public BasicDataReader<std::pair<size_t, IndigoFingerprint>, FingerprintTableReader, std::ifstream> {

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

        ReadValue readOne() override {
            uint64_t id;
            _binaryReader->read((char *) &id, sizeof id);
            IndigoFingerprint fingerprint;
            fingerprint.load(*_binaryReader);
            return {id, fingerprint};
        }

        bool eof() const override {
            return _fingerprintsInStream == 0;
        }
    };

} // qtr