#pragma once

#include "glog/logging.h"

#include "RawBucketIOConsts.h"
#include "basic_io/BasicDataWriter.h"

namespace qtr {

    // TODO class is not tested after refactoring
    class RawBucketWriter : public BasicDataWriter<raw_bucket_value_t, RawBucketWriter, std::ofstream> {
    private:
        uint64_t _writtenNumber;

    public:
        explicit RawBucketWriter(const std::filesystem::path &fileName) : BaseWriter(fileName), _writtenNumber(0) {
            LOG(INFO) << "Create raw bucket writer to " << fileName << " (" << _binaryWriter << ")";
            _binaryWriter->write((char *) &_writtenNumber, sizeof _writtenNumber); // reserve space for bucket size
        }

        ~RawBucketWriter() override {
            _binaryWriter->seekp(0, std::ios::beg);
            _binaryWriter->write((char *) &_writtenNumber, sizeof _writtenNumber); // write bucket size
            LOG(INFO) << "Delete raw bucket writer with " << _writtenNumber << " molecules (" << _binaryWriter << ")";
        }

        void write(const raw_bucket_value_t &value) override {
            _writtenNumber++;
            auto &[smiles, fingerprint] = value;
            fingerprint.saveBytes(*_binaryWriter);
            *_binaryWriter << smiles << '\n';
        }

        using BaseWriter::write;
    };

} // namespace qtr
