#pragma once

#include "glog/logging.h"

#include "basic_io/BasicDataWriter.h"
#include "fingerprint_table_io/FingerprintTableIOConsts.h"
#include "Fingerprint.h"

namespace qtr {

    class FingerprintTableWriter
            : public BasicDataWriter<std::pair<uint64_t, Fingerprint>, FingerprintTableWriter, std::ofstream> {
    private:
        uint64_t _writtenFingerprints;
        uint64_t _fingerprintSize;

    public:
        explicit FingerprintTableWriter(const std::filesystem::path &fileName) : BaseWriter(fileName),
                                                                                 _writtenFingerprints(0),
                                                                                 _fingerprintSize(0) {
            _binaryWriter->write((char *) &_writtenFingerprints,
                                 sizeof _writtenFingerprints); // reserve space for table size
            _binaryWriter->write((char *) &_fingerprintSize,
                                 sizeof _fingerprintSize); // reserve space for fingerprint size

            LOG(INFO) << "Create fingerprint table writer to " << fileName << " (" << _binaryWriter << ")";
        }

        ~FingerprintTableWriter() override {
            _binaryWriter->seekp(0, std::ios::beg);
            _binaryWriter->write((char *) &_writtenFingerprints, sizeof _writtenFingerprints);
            _binaryWriter->write((char *) &_fingerprintSize, sizeof _fingerprintSize);
            LOG(INFO) << "Delete fingerprint table writer with " << _writtenFingerprints << " molecules ("
                      << _binaryWriter << ")";
        }

        FingerprintTableWriter &operator<<(const WriteValue &value) override {
            auto &[id, fingerprint] = value;
            if (_writtenFingerprints == 0)
                _fingerprintSize = fingerprint.size();
            else
                assert(_fingerprintSize == fingerprint.size());
            _binaryWriter->write((char *) &id, sizeof id);
            fingerprint.dump(*_binaryWriter);
            _writtenFingerprints++;
            return *this;
        }

        using BaseWriter::operator<<;
    };

} // qtr