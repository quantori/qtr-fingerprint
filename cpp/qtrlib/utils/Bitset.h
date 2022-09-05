#pragma once

#include <cstddef>
#include <type_traits>
#include <fstream>
#include <iostream>
#include <cstdlib>

#include "Utils.h"

namespace qtr {

    // TODO test this class
    template<typename T, size_t S>
    class Bitset {
        static_assert(std::is_unsigned_v<T>());
        static_assert(sizeof(T) >= sizeof(char));

    protected:
        static const size_t _type_bits = fromBytesToBits(sizeof(T));
        static const size_t _data_length = divideIntegersCeil(S, _type_bits);
        static const size_t _size_in_bytes = divideIntegersCeil(S, BIT_IN_BYTE);
        T _data[_data_length];
    public:
        template<typename BinaryWriter>
        void dump(BinaryWriter &writer) {
            writer.write((char *) _data, _size_in_bytes);
        }

        template<typename BinaryReader>
        void load(BinaryReader &reader) {
            _data[_data_length - 1] = 0; // init extra bits with zeros
            reader.read(_data, _size_in_bytes);
        }

        bool operator<=(const Bitset &other) const {
            bool answer = true;
            for (size_t i = 0; i < _data_length && answer; i++) {
                answer &= (_data[i] & other._data[i]) == _data[i];
            }
            return answer;
        }

        bool get(size_t i) const {
            return _data[i >> _type_bits] >> lowerOrderBits(i, _type_bits) & T(1);
        }

        void set(size_t i, bool val) {
            T &num = _data[i >> _type_bits];
            size_t pos_in_num = lowerOrderBits(i, _type_bits);
            num &= T(-1) ^ (T(1) << pos_in_num);
            num |= T(val) << pos_in_num;
        }

    };

} // qrt