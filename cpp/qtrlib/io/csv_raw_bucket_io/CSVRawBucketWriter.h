#pragma once

#include "CSVRawBucketIOConsts.h"
#include "basic_io/BasicWriter.h"

namespace qtr {

    // TODO class is not tested after refactoring
    class CSVRawBucketWriter : public BasicWriter<csv_raw_bucket_value_t, CSVRawBucketWriter> {
    public:
        explicit CSVRawBucketWriter(std::ostream *stream) : BaseWriter(stream) {
            *_outStream << csvRawBucketHeader();
            LOG(INFO) << "Create CSV raw bucket writer (" << _stream << ")";
        };

        explicit CSVRawBucketWriter(const std::filesystem::path &fileName)
                : CSVRawBucketWriter(new std::ofstream(fileName)) {}

        CSVRawBucketWriter(const CSVRawBucketWriter &bucketWriter) = delete;

        ~CSVRawBucketWriter() override {
            LOG(INFO) << "Delete CSV raw bucket writer with " << _writtenNumber << " molecules (" << _stream << ")";
        }

        void write(const csv_raw_bucket_value_t &value) override {
            _writtenNumber++;
            const auto &[smiles, fingerprint] = value;
            *_outStream << smiles << csvSplitSymbol;
            for (size_t i = 0; i < IndigoFingerprint::sizeInBits; i++) {
                *_outStream << fingerprint[i];
                *_outStream << (i + 1 == IndigoFingerprint::sizeInBits ? '\n' : csvSplitSymbol);
            }
        }

        using BaseWriter::write;

    private:
        std::ostream *_outStream;
        uint64_t _writtenNumber;
    };

} // namespace qtr
