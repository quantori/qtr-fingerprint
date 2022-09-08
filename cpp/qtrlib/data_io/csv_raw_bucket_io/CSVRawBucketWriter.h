#pragma once

#include "CSVRawBucketIOConsts.h"
#include "basic_io/BasicDataWriter.h"

namespace qtr {

    // TODO class is not tested after refactoring
    class CSVRawBucketWriter : public BasicDataWriter<csv_raw_bucket_value_t, CSVRawBucketWriter, std::ofstream> {
    public:

        explicit CSVRawBucketWriter(const std::filesystem::path &fileName) : BaseWriter(fileName), _writtenNumber(0) {
            *_binaryWriter << csvRawBucketHeader();
            LOG(INFO) << "Create CSV raw bucket writer (" << _binaryWriter << ")";
        }

        CSVRawBucketWriter(const CSVRawBucketWriter &bucketWriter) = delete;

        ~CSVRawBucketWriter() override {
            LOG(INFO) << "Delete CSV raw bucket writer with " << _writtenNumber << " molecules (" << _binaryWriter
                      << ")";
        }

        void write(const csv_raw_bucket_value_t &value) override {
            _writtenNumber++;
            const auto &[smiles, fingerprint] = value;
            *_binaryWriter << smiles << csvSplitSymbol;
            for (size_t i = 0; i < IndigoFingerprint::sizeInBits; i++) {
                *_binaryWriter << fingerprint[i];
                *_binaryWriter << (i + 1 == IndigoFingerprint::sizeInBits ? '\n' : csvSplitSymbol);
            }
        }

        using BaseWriter::write;

    private:
        uint64_t _writtenNumber;
    };

} // namespace qtr
