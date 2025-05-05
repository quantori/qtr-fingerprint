#pragma once

#include <vector>
#include <cassert>
#include <stdexcept>
#include <limits>
#include <type_traits>

template<typename T>
class Bitset {
    static_assert(std::is_unsigned_v<T>, "T must be an unsigned type");
public:
    using BlockT = T;
protected:
    static constexpr size_t bitsPerBlock = std::numeric_limits<BlockT>::digits;

    std::vector<BlockT> _data;
    size_t _size;

    class Proxy {
    private:
        size_t _position;
        BlockT &_storage;

        static constexpr BlockT _makePosMask(size_t pos) {
            return BlockT(1) << pos;
        }

        void _setValue(bool value) {
            const BlockT mask = _makePosMask(_position);
            if (value) {
                _storage |= mask;
            } else {
                _storage &= ~mask;
            }
        }

    public:
        Proxy(size_t position, BlockT &storage) : _position(position), _storage(storage) {}

        operator bool() const {
            return (_storage & _makePosMask(_position)) != 0;
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

    static size_t bitSizeToBlockSize(size_t bitSize) {
        return (bitSize + bitsPerBlock - 1) / bitsPerBlock;
    }

public:
    explicit Bitset(size_t size)
            : _data(bitSizeToBlockSize(size), BlockT(0)),
              _size(size) {}

    Bitset() : _data(), _size(0) {}

    Bitset(const Bitset &other) = default;

    Bitset(Bitset &&other) noexcept = default;

    Bitset &operator=(const Bitset &other) = default;

    Bitset &operator=(Bitset &&other) noexcept = default;

    ~Bitset() = default;

    [[nodiscard]] size_t size() const {
        return _size;
    }

    BlockT *data() {
        return _data.data();
    }

    const BlockT *data() const {
        return _data.data();
    }

    [[nodiscard]] size_t dataVectorSize() const {
        return _data.size();
    }

    Bitset operator|(const Bitset &other) const {
        assert(_size == other._size && "Bitset sizes must match for operator|");
        Bitset answer(_size);
        for (size_t i = 0; i < _data.size(); ++i) {
            answer._data[i] = _data[i] | other._data[i];
        }
        return answer;
    }

    Bitset &operator|=(const Bitset &other) {
        assert(_size == other._size && "Bitset sizes must match for operator|=");
        for (size_t i = 0; i < _data.size(); ++i) {
            _data[i] |= other._data[i];
        }
        return *this;
    }

    bool operator==(const Bitset &other) const {
        if (_size != other._size) {
            return false;
        }
        return _data == other._data;
    }

    bool operator!=(const Bitset &other) const {
        if (_size != other._size) {
            return true;
        }
        return _data != other._data;
    }


    bool operator[](size_t i) const {
        assert(i < _size && "Bitset index out of range");

        size_t blockIdx = i / bitsPerBlock;
        size_t bitIdx = i % bitsPerBlock;

        return _data[blockIdx] & (BlockT(1) << bitIdx);
    }

    Proxy operator[](size_t i) {
        assert(i < _size && "Bitset index out of range");

        size_t blockIdx = i / bitsPerBlock;
        size_t bitIdx = i % bitsPerBlock;

        return Proxy(bitIdx, _data[blockIdx]);
    }
};