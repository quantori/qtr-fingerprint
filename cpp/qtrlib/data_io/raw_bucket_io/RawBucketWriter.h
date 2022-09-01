#pragma once

#include "glog/logging.h"

#include "RawBucketIOConsts.h"
#include "basic_io/BasicDataWriter.h"

namespace qtr {

    // TODO class is not tested after refactoring
    class RawBucketWriter : public BasicDataWriter<raw_bucket_value_t, RawBucketWriter> {
    private:
        uint64_t _writtenNumber;

    public:
        explicit RawBucketWriter(std::ostream *stream) : BaseWriter(stream), _writtenNumber(0) {
            _stream->write((char *) &_writtenNumber, sizeof _writtenNumber); // reserve space for bucket size
        };

        explicit RawBucketWriter(const std::filesystem::path &fileName) : RawBucketWriter(new std::ofstream(fileName)) {
            LOG(INFO) << "Create raw bucket writer to " << fileName << " (" << _stream << ")";
        }

        ~RawBucketWriter() override {
            _stream->seekp(0, std::ios::beg);
            _stream->write((char *) &_writtenNumber, sizeof _writtenNumber); // write bucket size
            LOG(INFO) << "Delete raw bucket writer with " << _writtenNumber << " molecules (" << _stream << ")";
        }

        void write(const raw_bucket_value_t &value) override {
            _writtenNumber++;
            auto &[smiles, fingerprint] = value;
            fingerprint.saveBytes(*_stream);
            *_stream << smiles << '\n';
        }

        using BaseWriter::write;
    };

} // namespace qtr
