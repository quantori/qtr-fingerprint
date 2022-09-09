#pragma once

#include <cstddef>
#include <type_traits>
#include <fstream>
#include <iostream>
#include <cstdlib>

#include "Utils.h"
#include "Fingerprint.h"

namespace qtr {

    using standardBitsetDataType = unsigned long long;

    // TODO test this class
    template<size_t S, typename T = standardBitsetDataType>
    class Bitset {
        static_assert(std::is_unsigned_v<T>);
        static_assert(sizeof(T) >= sizeof(char));

    protected:
        static const size_t _type_bits = fromBytesToBits(sizeof(T));
        static const size_t _data_length = divideIntegersCeil(S, _type_bits);
        static const size_t _index_shift = __builtin_clz(_type_bits);
        static const size_t sizeInBytes = divideIntegersCeil(S, BIT_IN_BYTE);
        T _data[_data_length];

        class Proxy {
        private:
            size_t _position;
            T &_storage;

            void _setValue(bool value) {
                _storage &= T(-1) ^ (T(1) << _position);
                _storage |= T(value) << _position;
            }

        public:
            Proxy(size_t position, T &storage) : _position(position), _storage(storage) {}

            operator bool() const {
                return _storage & (T(1) << _position);
            }

            Proxy &operator=(bool value) {
                _setValue(value);
                return *this;
            }

            Proxy &operator=(const Proxy &other) {
                _setValue(bool(other));
                return *this;
            }
        };

    public:
        Bitset() = default;

        Bitset(const Bitset &other) = default;

        static constexpr size_t size() {
            return S;
        }

        template<typename BinaryWriter>
        void dump(BinaryWriter &writer) const {
            writer.write((char *) _data, sizeInBytes);
        }

        template<typename BinaryReader>
        void load(BinaryReader &reader) {
            _data[_data_length - 1] = 0; // init extra bits with zeros
            reader.read((char *) _data, sizeInBytes);
        }

        bool operator<=(const Bitset &other) const {
            bool answer = true;
            for (size_t i = 0; i < _data_length && answer; i++) {
                answer &= (_data[i] & other._data[i]) == _data[i];
            }
            return answer;
        }

        // todo: test this function
        Bitset &reset() {
            std::memset(_data, 0, sizeof _data);
            return *this;
        }

        // todo: test this function
        Bitset operator|(const Bitset &other) const {
            Bitset answer;
            for (size_t i = 0; i < _data_length; i++) {
                answer._data[i] = _data[i] | other._data[i];
            }
        }

        bool operator==(const Bitset &other) const {
            bool answer = true;
            for (size_t i = 0; i < _data_length && answer; i++) {
                answer &= (_data[i] == other._data[i]);
            }
            return answer;
        }

        // todo: test this function
        Bitset &operator|=(const Bitset &other) {
            for (size_t i = 0; i < _data_length; i++) {
                _data[i] |= other._data[i];
            }
            return *this;
        }

        bool operator[](size_t i) const {
            return _data[i >> _index_shift] >> lowerOrderBits(i, _index_shift) & T(1);
        }

        Bitset::Proxy operator[](size_t i) {
            return Proxy(lowerOrderBits(i, _index_shift), _data[i >> _index_shift]);
        }
    };

} // qrt