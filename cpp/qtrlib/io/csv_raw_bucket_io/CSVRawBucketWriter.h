#pragma once

#include <iostream>
#include <filesystem>
#include <fstream>

#include "CSVRawBucketIO.h"

namespace qtr {

    class CSVRawBucketWriter {
    public:
        explicit CSVRawBucketWriter(std::ostream *outStream);

        explicit CSVRawBucketWriter(const std::filesystem::path &fileName);

        CSVRawBucketWriter(const CSVRawBucketWriter &bucketWriter) = delete;

        ~CSVRawBucketWriter();

        class Iterator {
        private:
            struct Proxy {
                Proxy(CSVRawBucketWriter::Iterator &it) : _iterator(it) {};

                Proxy &operator=(const csv_raw_bucket_value_t &value);

                CSVRawBucketWriter::Iterator &_iterator;
            };

            CSVRawBucketWriter *_writer;
            bool _isWritten;

        public:
            using iterator_category = std::output_iterator_tag;
            using difference_type = void;
            using value_type = csv_raw_bucket_value_t;
            using pointer = value_type *;
            using reference = value_type &;

            Iterator() : _writer(nullptr), _isWritten(false) {};

            explicit Iterator(CSVRawBucketWriter *writer);

            Iterator(const Iterator &it) = default;

            ~Iterator() = default;

            Proxy operator*();

            Iterator operator++();

            bool operator!=(const Iterator &it) const;

            bool isEnd() const;
        };

        void write(const csv_raw_bucket_value_t &value);

        void write(const std::vector<csv_raw_bucket_value_t> &values);

        Iterator begin();

        static Iterator end();

    private:
        std::ostream *_outStream;
        uint64_t _writtenNumber;
    };

} // namespace qtr
