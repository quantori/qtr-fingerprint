#include "CSVRawBucketReader.h"

#include "cassert"

namespace qtr {


    CSVRawBucketReader::CSVRawBucketReader(std::istream *inStream) {
        int c;
        auto header = csvRawBucketHeader();
        for (size_t i = 0; i < csvRawBucketHeaderSize; i++) {
            c = inStream->get();
            assert(c != EOF && c == header[i]);
        }
    }

    CSVRawBucketReader::CSVRawBucketReader(const std::filesystem::path &fileName) : CSVRawBucketReader(
            new std::ifstream(fileName)) {
        LOG(INFO) << "Create csv raw bucket reader from " << fileName << " (" << _inStream << ")";
    }

    CSVRawBucketReader::~CSVRawBucketReader() {
        LOG(INFO) << "Delete csv raw bucket reader (" << _inStream << ")";
        delete _inStream;
    }

    csv_raw_bucket_value_t CSVRawBucketReader::readOne() {
        int c;
        std::string smiles;
        while ((c = _inStream->get()) != csvSplitSymbol) {
            assert(c != EOF);
            smiles += (char)c;
        }
        IndigoFingerprint fingerprint;
        for (size_t i = 0; i < IndigoFingerprint::sizeInBits; i++) {
            c = _inStream->get();
            fingerprint[i] = bool(c - '0');
            c = _inStream->get();
            assert(c == csvSplitSymbol || i + 1 == IndigoFingerprint::sizeInBits);
            assert(c == '\n' || i + 1 != IndigoFingerprint::sizeInBits);
        }
        return {smiles, fingerprint};
    }

    std::vector<csv_raw_bucket_value_t> CSVRawBucketReader::readAll() {
        std::vector<csv_raw_bucket_value_t> csvRawBucket;
        std::copy(begin(), end(), std::back_inserter(csvRawBucket));
        return csvRawBucket;
    }

    CSVRawBucketReader::Iterator CSVRawBucketReader::begin() {
        return Iterator(this);
    }

    CSVRawBucketReader::Iterator CSVRawBucketReader::end() {
        return {};
    }

    CSVRawBucketReader::Iterator::Iterator(CSVRawBucketReader *reader) : _reader(reader), _isRead(false) {
        assert(reader != nullptr);
    }

    bool CSVRawBucketReader::Iterator::isEnd() const {
        return _reader == nullptr || _reader->_inStream->eof();
    }

    bool CSVRawBucketReader::Iterator::operator!=(const CSVRawBucketReader::Iterator &it) const {
        return (isEnd() ^ it.isEnd()) || (!isEnd() && !it.isEnd() && it._reader != _reader);
    }

    CSVRawBucketReader::Iterator::value_type CSVRawBucketReader::Iterator::operator*() {
        assert(_reader != nullptr && !_isRead && "Incorrect read of iterator");
        _isRead = true;
        return _reader->readOne();
    }

    CSVRawBucketReader::Iterator CSVRawBucketReader::Iterator::operator++() {
        _isRead = false;
        return *this;
    }
} // namespace qtr