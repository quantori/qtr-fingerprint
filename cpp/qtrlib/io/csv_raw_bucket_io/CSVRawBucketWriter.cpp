#include "CSVRawBucketWriter.h"

namespace qtr {


    CSVRawBucketWriter::CSVRawBucketWriter(std::ostream *outStream) : _outStream(outStream), _writtenNumber(0) {
        *_outStream << "SMILES" << csvSplitSymbol;
        for (size_t i = 0; i < IndigoFingerprint::sizeInBits; i++) {
            *_outStream << i;
            *_outStream << (i + 1 == IndigoFingerprint::sizeInBits ? '\n' : csvSplitSymbol);
        }
    }

    CSVRawBucketWriter::CSVRawBucketWriter(const std::filesystem::path &fileName) :
            CSVRawBucketWriter(new std::ofstream(fileName)) {}

    CSVRawBucketWriter::~CSVRawBucketWriter() {
        delete _outStream;
    }

    void CSVRawBucketWriter::write(const csv_raw_bucket_value_t &value) {
        _writtenNumber++;
        const auto &[smiles, fingerprint] = value;
        *_outStream << smiles << csvSplitSymbol;
        for (size_t i = 0; i < IndigoFingerprint::sizeInBits; i++) {
            *_outStream << fingerprint[i];
            *_outStream << (i + 1 == IndigoFingerprint::sizeInBits ? '\n' : csvSplitSymbol);
        }
    }

    void CSVRawBucketWriter::write(const std::vector<csv_raw_bucket_value_t> &values) {
        for (const auto &value: values) {
            write(value);
        }
    }


} // namespace qtr