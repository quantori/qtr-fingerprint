#pragma once

#include <iostream>
#include <fstream>
#include <filesystem>
#include <string>
#include <cassert>
#include <type_traits>

#include <glog/logging.h>


namespace qtr {

    template<typename Reader>
    class ReaderIterator {
    public:
        using iterator_category = std::input_iterator_tag;
        using difference_type = void;
        using value_type = typename Reader::ReadValue;
        using pointer = value_type *;
        using reference = value_type &;

        ReaderIterator() : _reader(nullptr), _isRead(false) {}

        explicit ReaderIterator(Reader *reader) : _reader(reader), _isRead(false) {
            assert(reader != nullptr);
        }

        ReaderIterator(const ReaderIterator &it) = default;

        ~ReaderIterator() = default;

        bool operator!=(const ReaderIterator &it) const {
            return (isEof() ^ it.isEof()) || (!isEof() && !it.isEof() && _reader != it._reader);
        }

        bool isEof() const {
            return _reader == nullptr || _reader->_binaryReader->peek() == EOF;
        }

        value_type operator*() {
            assert(_reader != nullptr && !_isRead && "Incorrect read of iterator");
            _isRead = true;
            value_type result;
            *_reader >> result;
            return result;
        }

        ReaderIterator &operator++() {
            _isRead = false;
            return *this;
        }

    protected:
        Reader *_reader;
        bool _isRead;
    };

    template<typename T, typename DataReader, typename BinaryReader>
    class BasicDataReader {
    public:
        using ReadValue = T;
        using BaseReader = BasicDataReader<T, DataReader, BinaryReader>;

        explicit BasicDataReader(const std::filesystem::path &filePath)
                : _binaryReader(new BinaryReader(filePath.c_str())) {}

        BasicDataReader(const BasicDataReader &basicReader) = delete;

        virtual ~BasicDataReader() {
            delete _binaryReader;
        }

        virtual ReaderIterator<DataReader> begin() {
            return ReaderIterator(dynamic_cast<DataReader *>(this));
        }

        virtual ReaderIterator<DataReader> end() {
            return {};
        }

        virtual BasicDataReader &operator>>(ReadValue &readValue) = 0;

        virtual BasicDataReader &operator>>(std::vector<ReadValue> &readValues) {
            std::copy(begin(), end(), std::back_inserter(readValues));
            return *this;
        }

        friend class ReaderIterator<DataReader>;

    protected:
        BinaryReader *_binaryReader;
    };

} // qtr

