#pragma once

#include <string>

#include "SmilesTableIOConsts.h"
#include "basic_io/BasicDataReader.h"

namespace qtr {

    class SmilesTableReader
            : public BasicDataReader<smiles_table_value_t, SmilesTableReader, std::ifstream> {
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

        SmilesTableReader &operator>>(ReadValue &value) override {
            auto &[id, smiles] = value;
            _binaryReader->read((char *) &id, sizeof id);
            char symbol;
            while ((symbol = (char) _binaryReader->get()) != '\n') {
                smiles += symbol;
            }
            return *this;
        }

        using BaseReader::operator>>;
    };


} // qtr