#pragma once

#include "RawBucketIOConsts.h"
#include "basic_io/BasicDataReader.h"

namespace qtr {

    // TODO class is not tested after refactoring
    class RawBucketReader : public BasicDataReader<raw_bucket_value_t, RawBucketReader, std::ifstream> {
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

        raw_bucket_value_t readOne() override {
            assert(_moleculesInStream != 0);
            _moleculesInStream--;
            IndigoFingerprint fingerprint;
            std::string smiles;
            fingerprint.load(*_binaryReader);
            char symbol;
            while ((symbol = (char) _binaryReader->get()) != '\n') {
                smiles += symbol;
            }
            return {smiles, fingerprint};
        }

        bool eof() const override {
            return _moleculesInStream == 0;
        }
    };

}  // namespace qtr
