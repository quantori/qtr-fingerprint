#pragma once

#include <string>

#include "basic_io/BasicDataReader.h"

namespace qtr {

    // TODO: test this class
    class SmilesTableReader : public BasicDataReader<std::pair<uint64_t, std::string>, SmilesTableReader> {
    private:
        uint64_t _smilesInStream;

    public:
        explicit SmilesTableReader(std::istream *stream) : BaseReader(stream), _smilesInStream(0) {
            _stream->read((char *) &_smilesInStream, sizeof _smilesInStream);
        }

        explicit SmilesTableReader(const std::filesystem::path &fileName)
                : SmilesTableReader(new std::ifstream(fileName)) {
            LOG(INFO) << "Create SMILES table reader with " << _smilesInStream << " SMILES (" << _stream << ")";
        }

        ~SmilesTableReader() override {
            LOG(INFO) << "Delete SMILES table reader (" << _stream << ")";
        }

        ReadValue readOne() override {
            uint64_t id;
            _stream->read((char *) &id, sizeof id);
            std::string smiles;
            char symbol;
            while ((symbol = (char) _stream->get()) != '\n') {
                smiles += symbol;
            }
            return {id, smiles};
        }

        bool isEof() const override {
            return _smilesInStream == 0;
        }
    };


} // qtr