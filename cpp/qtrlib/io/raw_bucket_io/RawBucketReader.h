#pragma once

#include "RawBucketIOConsts.h"
#include "basic_io/BasicReader.h"

namespace qtr {

    // TODO class is not tested after refactoring
    class RawBucketReader : public BasicReader<raw_bucket_value_t, RawBucketReader> {
    private:
        uint64_t _moleculesInStream;

    public:
        explicit RawBucketReader(std::istream *stream) : BaseReader(stream), _moleculesInStream(0) {
            _stream->read((char *) &_moleculesInStream, sizeof _moleculesInStream);
        }

        explicit RawBucketReader(const std::filesystem::path &fileName) : RawBucketReader(new std::ifstream(fileName)) {
            LOG(INFO) << "Create raw bucket reader from " << fileName << " with " << _moleculesInStream
                      << " molecules (" << _stream << ")";
        }

        ~RawBucketReader() override {
            LOG(INFO) << "Delete raw bucket reader (" << _stream << ")";
        }

        raw_bucket_value_t readOne() override {
            assert(_moleculesInStream != 0);
            _moleculesInStream--;
            IndigoFingerprint fingerprint;
            std::string smiles;
            fingerprint.readFrom(*_stream);
            char symbol;
            while ((symbol = (char) _stream->get()) != '\n') {
                smiles += symbol;
            }
            return {smiles, fingerprint};
        }

        bool isEof() const override {
            return _moleculesInStream == 0;
        }
    };

}  // namespace qtr
