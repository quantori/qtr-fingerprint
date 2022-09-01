#pragma once

#include "CSVRawBucketIOConsts.h"
#include "basic_io/BasicDataReader.h"

namespace qtr {

    // TODO class is not tested after refactoring
    class CSVRawBucketReader : public BasicDataReader<csv_raw_bucket_value_t, CSVRawBucketReader> {
    public:
        explicit CSVRawBucketReader(std::istream *stream): BaseReader(stream) {
            int c;
            auto header = csvRawBucketHeader();
            for (size_t i = 0; i < csvRawBucketHeaderSize; i++) {
                c = _stream->get();
                assert(c != EOF && c == header[i]);
            }
        };

        explicit CSVRawBucketReader(const std::filesystem::path &fileName)
                : CSVRawBucketReader(new std::ifstream(fileName)) {
            LOG(INFO) << "Create csv raw bucket reader from " << fileName << " (" << _stream << ")";
        }

        CSVRawBucketReader(const CSVRawBucketReader &bucketLoader) = delete;

        ~CSVRawBucketReader() override {
            LOG(INFO) << "Delete csv raw bucket reader (" << _stream << ")";
        }

        csv_raw_bucket_value_t readOne() override {
            int c;
            std::string smiles;
            while ((c = _stream->get()) != csvSplitSymbol) {
                assert(c != EOF);
                smiles += (char) c;
            }
            IndigoFingerprint fingerprint;
            for (size_t i = 0; i < IndigoFingerprint::sizeInBits; i++) {
                c = _stream->get();
                fingerprint[i] = bool(c - '0');
                c = _stream->get();
                assert(c == csvSplitSymbol || i + 1 == IndigoFingerprint::sizeInBits);
                assert(c == '\n' || i + 1 != IndigoFingerprint::sizeInBits);
            }
            return {smiles, fingerprint};
        }

        bool isEof() const override {
            return _stream->eof();
        }
    };

} // namespace qtr