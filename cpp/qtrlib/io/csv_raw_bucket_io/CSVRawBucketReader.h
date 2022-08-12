#pragma once

#include <fstream>
#include <filesystem>

#include "CSVRawBucketIOConsts.h"

namespace qtr {

    class CSVRawBucketReader {
    public:
        explicit CSVRawBucketReader(std::istream *inStream);

        explicit CSVRawBucketReader(const std::filesystem::path &fileName);

        CSVRawBucketReader(const CSVRawBucketReader &bucketLoader) = delete;

        ~CSVRawBucketReader();

        class Iterator {
        public:
            using iterator_category = std::input_iterator_tag;
            using difference_type = void;
            using value_type = qtr::csv_raw_bucket_value_t;
            using pointer = value_type *;
            using reference = value_type &;

            Iterator() : _reader(nullptr), _isRead(false) {};

            explicit Iterator(CSVRawBucketReader *reader);

            Iterator(const Iterator &it) = default;

            ~Iterator() = default;

            bool isEnd() const;

            bool operator!=(const Iterator &it) const;

            value_type operator*();

            Iterator operator++();

        private:
            CSVRawBucketReader *_reader;
            bool _isRead;
        };

        csv_raw_bucket_value_t readOne();

        std::vector<csv_raw_bucket_value_t> readAll();

        Iterator begin();

        static Iterator end();

    private:
        uint64_t _moleculesInStream;
        std::istream *_inStream;
    };

} // namespace qtr