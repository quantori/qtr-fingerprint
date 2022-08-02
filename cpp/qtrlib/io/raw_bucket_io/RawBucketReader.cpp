#include "RawBucketReader.h"

namespace qtr {

    RawBucketReader::Iterator RawBucketReader::begin() {
        return Iterator(this);
    }


    RawBucketReader::Iterator RawBucketReader::end() {
        return {};
    }

    RawBucketReader::RawBucketReader(std::istream *inStream) : _inStream(inStream), _bucketsInStream(0) {
        _inStream->read((char *) &_bucketsInStream, sizeof _bucketsInStream);
    }

    RawBucketReader::Iterator::value_type RawBucketReader::readOne() {
        assert(_bucketsInStream != 0);
        _bucketsInStream--;
        IndigoFingerprint fingerprint;
        std::string smiles;
        fingerprint.readFrom(*_inStream);
        char symbol;
        while ((symbol = (char) _inStream->get()) != '\n') {
            smiles += symbol;
        }
        return {smiles, fingerprint};
    }

    RawBucketReader::~RawBucketReader() {
        delete _inStream;
    }

    std::vector<raw_bucket_value_t> RawBucketReader::readAll() {
        std::vector<raw_bucket_value_t> rawBucket;
        std::copy(this->begin(), this->end(), std::back_inserter(rawBucket));
        return rawBucket;
    }

    bool RawBucketReader::Iterator::isEnd() const {
        return _reader == nullptr || _reader->_bucketsInStream == 0;
    }

    RawBucketReader::Iterator::Iterator(RawBucketReader *reader) : _reader(reader), _isRead(false) {
        assert(reader != nullptr);
    }

    bool RawBucketReader::Iterator::operator!=(const RawBucketReader::Iterator &it) const {
        return (it.isEnd() ^ isEnd()) || (!it.isEnd() && !isEnd() && it._reader != _reader);
    }

    RawBucketReader::Iterator::value_type RawBucketReader::Iterator::operator*() {
        assert(_reader != nullptr && !_isRead && "Incorrect read of iterator");
        _isRead = true;
        return _reader->readOne();
    }

    RawBucketReader::Iterator RawBucketReader::Iterator::operator++() {
        _isRead = false;
        return *this;
    }

}