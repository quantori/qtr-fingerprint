#pragma once

#include "RawBucketIOConsts.h"
#include "basic_io/BasicDataReader.h"
#include "io/BufferedReader.h"

namespace qtr {

    class RawBucketReader : public BasicDataReader<raw_bucket_value_t, RawBucketReader, BufferedReader<>> {
    private:
        uint64_t _moleculesInStream;

    public:

        explicit RawBucketReader(const std::filesystem::path &fileName) : BaseReader(fileName), _moleculesInStream(0) {
            _binaryReader->read((char *) &_moleculesInStream, sizeof _moleculesInStream);
            LOG(INFO) << "Create raw bucket reader from " << fileName << " with " << _moleculesInStream
                      << " molecules (" << _binaryReader << ")";
        }

        ~RawBucketReader() override {
            LOG(INFO) << "Delete raw bucket reader (" << _binaryReader << ")";
        }

        RawBucketReader &operator>>(raw_bucket_value_t &value) override {
            assert(_moleculesInStream != 0);
            _moleculesInStream--;
            auto &[smiles, fingerprint] = value;
            fingerprint.load(*_binaryReader, 0);
            char symbol;
            while ((symbol = (char) _binaryReader->get()) != '\n') {
                smiles += symbol;
            }
            return *this;
        }

        using BaseReader::operator>>;
    };

}  // namespace qtr
