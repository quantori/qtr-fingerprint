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
            return _reader == nullptr || _reader->isEof();
        }

        value_type operator*() {
            assert(_reader != nullptr && !_isRead && "Incorrect read of iterator");
            _isRead = true;
            return _reader->readOne();
        }

        ReaderIterator &operator++() {
            _isRead = false;
            return *this;
        }

    protected:
        Reader *_reader;
        bool _isRead;
    };

    template<typename T, typename Reader>
    class BasicReader {
    public:
        using ReadValue = T;
        using BaseReader = BasicReader<T, Reader>;

        explicit BasicReader(std::istream *stream) : _stream(stream) {};

        explicit BasicReader(const std::filesystem::path &filePath) : BasicReader(new std::ifstream(filePath)) {}

        BasicReader(const BasicReader &basicReader) = delete;

        virtual ~BasicReader() {
            delete _stream;
        }

        virtual ReaderIterator<Reader> begin() {
            return ReaderIterator(dynamic_cast<Reader *>(this));
        }

        virtual ReaderIterator<Reader> end() {
            return {};
        }

        virtual ReadValue readOne() = 0;

        virtual std::vector<ReadValue> readAll() {
            std::vector<ReadValue> result;
            std::copy(begin(), end(), std::back_inserter(result));
            return result;
        }

        virtual bool isEof() const = 0;


        friend class ReaderIterator<Reader>;

    protected:
        std::istream *_stream;
    };

} // qtr

