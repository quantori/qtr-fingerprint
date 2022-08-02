#include "RawBucketWriter.h"

namespace qtr {
    RawBucketWriter::RawBucketWriter(std::ostream *outStream) : _outStream(outStream), _writtenNumber(0) {
        _outStream->write((char *) &_writtenNumber, sizeof _writtenNumber); // reserve space for bucket size
    }

    RawBucketWriter::~RawBucketWriter() {
        _outStream->seekp(0, std::ios::beg);
        _outStream->write((char *) &_writtenNumber, sizeof _writtenNumber); // write bucket size
        delete _outStream;
    }

    void RawBucketWriter::write(const raw_bucket_value_t &value) {
        _writtenNumber++;
        auto &[smiles, fingerprint] = value;
        fingerprint.saveBytes(*_outStream);
        *_outStream << smiles << '\n';
    }

    RawBucketWriter::Iterator RawBucketWriter::begin() {
        return Iterator(this);
    }

    RawBucketWriter::Iterator RawBucketWriter::end() {
        return {};
    }

    void RawBucketWriter::write(const std::vector<raw_bucket_value_t> &values) {
        std::copy(values.begin(), values.end(), this->begin());
    }

    RawBucketWriter::Iterator::Proxy &RawBucketWriter::Iterator::Proxy::operator=(const raw_bucket_value_t &value) {
        _iterator._writer->write(value);
        return *this;
    }

    RawBucketWriter::Iterator::Proxy RawBucketWriter::Iterator::operator*() {
        assert(!_isWritten);
        _isWritten = true;
        return {*this};
    }

    RawBucketWriter::Iterator::Iterator(RawBucketWriter *writer) : _writer(writer), _isWritten(false) {
        assert(writer != nullptr);
    }

    RawBucketWriter::Iterator RawBucketWriter::Iterator::operator++() {
        assert(_isWritten && "increment not written iterator");
        _isWritten = false;
        return *this;
    }

    bool RawBucketWriter::Iterator::isEnd() const {
        return _writer == nullptr;
    }

    bool RawBucketWriter::Iterator::operator!=(const RawBucketWriter::Iterator &it) const {
        return (isEnd() ^ it.isEnd()) || (!isEnd() && !it.isEnd() && _writer == it._writer);
    }
} // namespace qtr
