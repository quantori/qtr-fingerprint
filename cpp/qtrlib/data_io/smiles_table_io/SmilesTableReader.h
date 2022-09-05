#pragma once

#include <string>

#include "basic_io/BasicDataReader.h"

namespace qtr {

    // TODO: test this class
    class SmilesTableReader
            : public BasicDataReader<std::pair<uint64_t, std::string>, SmilesTableReader, std::ifstream> {
    private:
        uint64_t _smilesInStream;

    public:
        explicit SmilesTableReader(const std::filesystem::path &fileName)
                : BaseReader(fileName), _smilesInStream(0) {
            _binaryReader->read((char *) &_smilesInStream, sizeof _smilesInStream);
            LOG(INFO) << "Create SMILES table reader with " << _smilesInStream << " SMILES (" << _binaryReader << ")";
        }

        ~SmilesTableReader() override {
            LOG(INFO) << "Delete SMILES table reader (" << _binaryReader << ")";
        }

        ReadValue readOne() override {
            uint64_t id;
            _binaryReader->read((char *) &id, sizeof id);
            std::string smiles;
            char symbol;
            while ((symbol = (char) _binaryReader->get()) != '\n') {
                smiles += symbol;
            }
            return {id, smiles};
        }

        bool eof() const override {
            return _smilesInStream == 0;
        }
    };


} // qtr