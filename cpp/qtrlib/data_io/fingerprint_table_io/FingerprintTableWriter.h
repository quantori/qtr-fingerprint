#pragma once

#include "glog/logging.h"

#include "basic_io/BasicDataWriter.h"
#include "fingerprint_table_io/FingerprintTableIOConsts.h"
#include "Fingerprint.h"

namespace qtr {

    class FingerprintTableWriter
            : public BasicDataWriter<std::pair<uint64_t, IndigoFingerprint>, FingerprintTableWriter, std::ofstream> {
    private:
        uint64_t _writtenFingerprints;

    public:
        explicit FingerprintTableWriter(const std::filesystem::path &fileName) : BaseWriter(fileName),
                                                                                 _writtenFingerprints(0) {
            _binaryWriter->write((char *) &_writtenFingerprints,
                                 sizeof _writtenFingerprints); // reserve space for table size
            LOG(INFO) << "Create fingerprint table writer to " << fileName << " (" << _binaryWriter << ")";
        }

        ~FingerprintTableWriter() override {
            _binaryWriter->seekp(0, std::ios::beg);
            _binaryWriter->write((char *) &_writtenFingerprints, sizeof _writtenFingerprints); // write bucket size
            LOG(INFO) << "Delete fingerprint table writer with " << _writtenFingerprints << " molecules ("
                      << _binaryWriter << ")";
        }

        FingerprintTableWriter &operator<<(const WriteValue &value) override {
            _writtenFingerprints++;
            auto &[id, fingerprint] = value;
            _binaryWriter->write((char *) &id, sizeof id);
            fingerprint.dump(*_binaryWriter);
            return *this;
        }

        using BaseWriter::operator<<;
    };

} // qtr