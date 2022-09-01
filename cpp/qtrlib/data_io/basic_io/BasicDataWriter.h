#pragma once

#include <iostream>
#include <fstream>
#include <filesystem>
#include <vector>
#include <cassert>

namespace qtr {

    template<typename Writer>
    class WriterIterator {
    private:
        struct Proxy;

        Writer *_writer;
        bool _isWritten;

    public:
        using iterator_category = std::output_iterator_tag;
        using difference_type = void;
        using value_type = typename Writer::WriteValue;
        using pointer = value_type *;
        using reference = value_type &;

        WriterIterator() : _writer(nullptr), _isWritten(false) {};

        explicit WriterIterator(Writer *writer) : _writer(writer), _isWritten(false) {
            assert(writer != nullptr);
        }

        WriterIterator(const WriterIterator &it) = default;

        virtual ~WriterIterator() = default;

        virtual bool operator!=(const WriterIterator &it) const {
            return _writer != it._writer;
        }

        virtual Proxy operator*() {
            assert(!_isWritten);
            _isWritten = true;
            return {*this};
        }

        virtual WriterIterator &operator++() {
            assert(_isWritten && "increment not written iterator");
            _isWritten = false;
            return *this;
        }

    private:
        struct Proxy {
            Proxy(WriterIterator &it) : _iterator(it) {};

            Proxy &operator=(const value_type &value) {
                _iterator._writer->write(value);
                return *this;
            }

            WriterIterator &_iterator;
        };
    };

    template<typename T, typename Writer>
    class BasicDataWriter {
    public:
        using WriteValue = T;
        using BaseWriter = BasicDataWriter<T, Writer>;

        explicit BasicDataWriter(std::ostream *stream) : _stream(stream), _writtenSmiles(0) {}

        BasicDataWriter(const BasicDataWriter &bucketWriter) = delete;

        virtual ~BasicDataWriter() {
            delete _stream;
        }

        virtual WriterIterator<Writer> begin() {
            return WriterIterator(dynamic_cast<Writer *>(this));
        }

        virtual WriterIterator<Writer> end() {
            return WriterIterator<Writer>();
        }

        virtual void write(const WriteValue &value) = 0;

        virtual void write(const std::vector<WriteValue> &values) {
            std::copy(values.begin(), values.end(), begin());
        }

    protected:
        std::ostream *_stream;
        uint64_t _writtenSmiles;
    };
} // qtr