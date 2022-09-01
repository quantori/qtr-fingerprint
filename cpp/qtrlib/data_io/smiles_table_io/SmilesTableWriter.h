#pragma once

#include <string>

#include "glog/logging.h"

#include "basic_io/BasicDataWriter.h"

namespace qtr {

    // TODO test this class
    class SmilesTableWriter : public BasicDataWriter<std::pair<uint64_t, std::string>, SmilesTableWriter> {
    private:
        uint64_t _writtenSmiles;

    public:
        explicit SmilesTableWriter(std::ostream *stream) : BaseWriter(stream), _writtenSmiles(0) {
            _stream->write((char *) &_writtenSmiles, sizeof _writtenSmiles); // reserve space for table size
        }

        explicit SmilesTableWriter(const std::filesystem::path &fileName)
                : SmilesTableWriter(new std::ofstream(fileName)) {
            LOG(INFO) << "Create SMILES talbe writer to " << fileName << " (" << _stream << ")";
        }

        ~SmilesTableWriter() override {
            _stream->seekp(0, std::ios::beg);
            _stream->write((char *) &_writtenSmiles, sizeof _writtenSmiles); // write bucket size
            LOG(INFO) << "Delete SMILES table writer with " << _writtenSmiles << " molecules (" << _stream << ")";
        }

        void write(const WriteValue &value) override {
            _writtenSmiles++;
            auto &[id, smiles] = value;
            _stream->write((char *) &id, sizeof id);
            *_stream << smiles << '\n';
        }
    };

} // qtr