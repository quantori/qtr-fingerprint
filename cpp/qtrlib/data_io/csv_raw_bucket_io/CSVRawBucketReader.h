#pragma once

#include "CSVRawBucketIOConsts.h"
#include "basic_io/BasicDataReader.h"

namespace qtr {

    // TODO class is not tested after refactoring
    class CSVRawBucketReader : public BasicDataReader<csv_raw_bucket_value_t, CSVRawBucketReader, std::ifstream> {
    public:

        explicit CSVRawBucketReader(const std::filesystem::path &fileName) : BaseReader(fileName) {
            LOG(INFO) << "Create csv raw bucket reader from " << fileName << " (" << _binaryReader << ")";

            int c;
            auto header = csvRawBucketHeader();
            for (size_t i = 0; i < csvRawBucketHeaderSize; i++) {
                c = _binaryReader->get();
                assert(c != EOF && c == header[i]);
            }
        }

        CSVRawBucketReader(const CSVRawBucketReader &bucketLoader) = delete;

        ~CSVRawBucketReader() override {
            LOG(INFO) << "Delete csv raw bucket reader (" << _binaryReader << ")";
        }

        csv_raw_bucket_value_t readOne() override {
            int c;
            std::string smiles;
            while ((c = _binaryReader->get()) != csvSplitSymbol) {
                assert(c != EOF);
                smiles += (char) c;
            }
            IndigoFingerprint fingerprint;
            for (size_t i = 0; i < IndigoFingerprint::sizeInBits; i++) {
                c = _binaryReader->get();
                fingerprint[i] = bool(c - '0');
                c = _binaryReader->get();
                assert(c == csvSplitSymbol || i + 1 == IndigoFingerprint::sizeInBits);
                assert(c == '\n' || i + 1 != IndigoFingerprint::sizeInBits);
            }
            return {smiles, fingerprint};
        }

        bool eof() const override {
            return _binaryReader->eof();
        }
    };

} // namespace qtr