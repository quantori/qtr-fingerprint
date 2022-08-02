#pragma once

#include <iterator>
#include <algorithm>
#include <string>
#include <fstream>
#include <filesystem>
#include <cassert>

#include "RawBucketsIO.h"

namespace qtr {

    class RawBucketReader {
    public:
        explicit RawBucketReader(std::istream *inStream);

        explicit RawBucketReader(const std::filesystem::path &fileName) :
                RawBucketReader(new std::ifstream(fileName)) {};

        RawBucketReader(const RawBucketReader &bucketLoader) = delete;

        ~RawBucketReader();

        class Iterator {
        public:
            using iterator_category = std::input_iterator_tag;
            using difference_type = void;
            using value_type = qtr::raw_bucket_value_t;
            using pointer = value_type *;
            using reference = value_type &;

            Iterator() : _reader(nullptr), _isRead(false) {};

            explicit Iterator(RawBucketReader *reader);

            Iterator(const Iterator& it) = default;

            ~Iterator() = default;

            bool isEnd() const;

            bool operator!=(const Iterator &it) const;

            value_type operator*();

            Iterator operator++();

        private:
            RawBucketReader *_reader;
            bool _isRead;
        };

        raw_bucket_value_t readOne();

        std::vector<raw_bucket_value_t> readAll();

        Iterator begin();

        static Iterator end();

    private:
        uint64_t _bucketsInStream;
        std::istream *_inStream;
    };

}  // namespace qtr
