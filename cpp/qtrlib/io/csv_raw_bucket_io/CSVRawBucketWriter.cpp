#include "CSVRawBucketWriter.h"

#include <cassert>

namespace qtr {


    CSVRawBucketWriter::CSVRawBucketWriter(std::ostream *outStream) : _outStream(outStream), _writtenNumber(0) {
        *_outStream << csvRawBucketHeader();
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

    CSVRawBucketWriter::Iterator CSVRawBucketWriter::begin() {
        return Iterator(this);
    }

    CSVRawBucketWriter::Iterator CSVRawBucketWriter::end() {
        return {};
    }

    CSVRawBucketWriter::Iterator::Proxy &
    CSVRawBucketWriter::Iterator::Proxy::operator=(const csv_raw_bucket_value_t &value) {
        _iterator._writer->write(value);
        return *this;
    }

    CSVRawBucketWriter::Iterator::Iterator(CSVRawBucketWriter *writer) : _writer(writer), _isWritten(false) {
        assert(writer != nullptr);
    }

    CSVRawBucketWriter::Iterator::Proxy CSVRawBucketWriter::Iterator::operator*() {
        assert(!_isWritten);
        _isWritten = true;
        return {*this};
    }

    CSVRawBucketWriter::Iterator CSVRawBucketWriter::Iterator::operator++() {
        assert(_isWritten && "increment not written iterator");
        _isWritten = false;
        return *this;
    }

    bool CSVRawBucketWriter::Iterator::operator!=(const CSVRawBucketWriter::Iterator &it) const {
        return (isEnd() ^ it.isEnd()) || (!isEnd() && !it.isEnd() && _writer == it._writer);
    }

    bool CSVRawBucketWriter::Iterator::isEnd() const {
        return _writer == nullptr;
    }

} // namespace qtr