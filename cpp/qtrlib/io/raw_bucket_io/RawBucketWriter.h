#pragma once

#include <iostream>
#include <filesystem>
#include <fstream>

#include "RawBucketsIO.h"

namespace qtr {

    class RawBucketWriter {
    public:
        explicit RawBucketWriter(std::ostream *outStream);

        explicit RawBucketWriter(const std::filesystem::path &fileName) :
                RawBucketWriter(new std::ofstream(fileName)) {};

        RawBucketWriter(const RawBucketWriter &bucketWriter) = delete;

        ~RawBucketWriter();

        class Iterator {
        private:
            struct Proxy {
                Proxy(RawBucketWriter::Iterator &it) : _iterator(it) {};

                Proxy &operator=(const raw_bucket_value_t &value);

                RawBucketWriter::Iterator &_iterator;
            };

            RawBucketWriter *_writer;
            bool _isWritten;

        public:
            using iterator_category = std::output_iterator_tag;
            using difference_type = void;
            using value_type = raw_bucket_value_t;
            using pointer = value_type *;
            using reference = value_type &;

            Iterator() : _writer(nullptr), _isWritten(false) {};

            explicit Iterator(RawBucketWriter *writer);

            Iterator(const Iterator &it) = default;

            ~Iterator() = default;

            Proxy operator*();

            Iterator operator++();

            bool operator!=(const Iterator &it) const;

            bool isEnd() const;
        };

        void write(const raw_bucket_value_t &value);

        void write(const std::vector<raw_bucket_value_t> &values);

        Iterator begin();

        static Iterator end();

    private:
        uint64_t _writtenNumber;
        std::ostream *_outStream;
    };

} // namespace qtr