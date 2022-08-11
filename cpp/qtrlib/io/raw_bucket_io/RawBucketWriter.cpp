#include "RawBucketWriter.h"

namespace qtr {
    RawBucketWriter::RawBucketWriter(std::ostream *outStream) : _outStream(outStream), _writtenNumber(0) {
        _outStream->write((char *) &_writtenNumber, sizeof _writtenNumber); // reserve space for bucket size
    }

    RawBucketWriter::RawBucketWriter(const std::filesystem::path &fileName) :
            RawBucketWriter(new std::ofstream(fileName)) {
        LOG(INFO) << "Create raw bucket writer to " << fileName << " (" << _outStream << ")";
    }

    RawBucketWriter::~RawBucketWriter() {
        _outStream->seekp(0, std::ios::beg);
        _outStream->write((char *) &_writtenNumber, sizeof _writtenNumber); // write bucket size
        LOG(INFO) << "Delete raw bucket writer with " << _writtenNumber << " molecules (" << _outStream << ")";
        delete _outStream;
    }

    void RawBucketWriter::write(const raw_bucket_value_t &value) {
        _writtenNumber++;
        auto &[smiles, fingerprint] = value;
        fingerprint.saveBytes(*_outStream);
        *_outStream << smiles << '\n';
    }

    void RawBucketWriter::write(const std::vector<raw_bucket_value_t> &values) {
        std::copy(values.begin(), values.end(), this->begin());
    }

    RawBucketWriter::Iterator::Iterator(RawBucketWriter *writer) : _writer(writer), _isWritten(false) {
        assert(writer != nullptr);
    }

    RawBucketWriter::Iterator RawBucketWriter::begin() {
        return Iterator(this);
    }

    RawBucketWriter::Iterator RawBucketWriter::end() {
        return {};
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

    RawBucketWriter::Iterator RawBucketWriter::Iterator::operator++() {
        assert(_isWritten && "increment not written iterator");
        _isWritten = false;
        return *this;
    }

    bool RawBucketWriter::Iterator::operator!=(const RawBucketWriter::Iterator &it) const {
        return (isEnd() ^ it.isEnd()) || (!isEnd() && !it.isEnd() && _writer == it._writer);
    }

    bool RawBucketWriter::Iterator::isEnd() const {
        return _writer == nullptr;
    }

} // namespace qtr
