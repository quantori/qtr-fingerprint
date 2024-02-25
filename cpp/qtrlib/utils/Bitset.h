#pragma once

#include <cstddef>
#include <type_traits>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cmath>

#include "Utils.h"

namespace qtr {

    using standardBitsetDataType = unsigned long long;

    template<typename T = standardBitsetDataType>
    class Bitset {
        static_assert(std::is_unsigned_v<T>);
        static_assert(sizeof(T) >= sizeof(char));

    protected:
        static const size_t TypeSizeInBits = fromBytesToBits(sizeof(T));
        static const size_t IndexShift = log2Floor(TypeSizeInBits);
        static_assert((1ull << IndexShift) == TypeSizeInBits);
        std::vector<T> _data;
        size_t _size;

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

        [[nodiscard]] size_t sizeInBytes() const {
            return divideIntegersCeil(_size, BIT_IN_BYTE);
        }

    public:
        explicit Bitset(size_t size) : _data(divideIntegersCeil(size, TypeSizeInBits), T(0)), _size(size) {
            reset();
        }

        Bitset(const Bitset &other) = default;

        Bitset() = default;

        [[nodiscard]] size_t size() const {
            return _size;
        }

        template<typename BinaryWriter>
        void dump(BinaryWriter &writer) const {
            writer.write((char *) _data.data(), sizeInBytes());
        }

        template<typename BinaryReader>
        void load(BinaryReader &reader, size_t size) {
            assert(size != 0);
            new(this) Bitset(size);
            _data.back() = 0; // init extra bits with zeros
            reader.read((char *) _data.data(), sizeInBytes());
        }

        bool operator<=(const Bitset &other) const {
            assert(_size == other._size);
            for (size_t i = 0; i < _data.size(); i++) {
                if (_data[i] & ~other._data[i])
                    return false;
            }
            return true;
        }

        Bitset &reset() {
            std::memset(_data.data(), 0, sizeInBytes());
            return *this;
        }

        Bitset &set() {
            std::memset(_data.data(), 255, sizeInBytes());
            _data.back() &= T(1) << (size() % TypeSizeInBits) - 1;
            return *this;
        }

        T *data() {
            return _data.data();
        }

        Bitset operator|(const Bitset &other) const {
            assert(_size == other._size);
            Bitset answer(_size);
            for (size_t i = 0; i < _data.size(); i++) {
                answer._data[i] = _data[i] | other._data[i];
            }
            return answer;
        }

        bool operator==(const Bitset &other) const {
            assert(_size == other._size);
            return _data == other._data;
        }

        Bitset &operator|=(const Bitset &other) {
            assert(_size == other._size);
            for (size_t i = 0; i < _data.size(); i++) {
                _data[i] |= other._data[i];
            }
            return *this;
        }

        bool operator[](size_t i) const {
            return _data[i >> IndexShift] >> lowerOrderBits(i, IndexShift) & T(1);
        }

        Bitset::Proxy operator[](size_t i) {
            return Proxy(lowerOrderBits(i, IndexShift), _data[i >> IndexShift]);
        }
    };

} // qrt