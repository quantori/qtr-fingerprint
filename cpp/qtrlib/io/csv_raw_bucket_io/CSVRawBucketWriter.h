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

        void write(const csv_raw_bucket_value_t &value);

        void write(const std::vector<csv_raw_bucket_value_t> &values);

    private:
        std::ostream *_outStream;
        uint64_t _writtenNumber;
    };

} // namespace qtr
