#pragma once

#include "MemoryMapBytesStorage.h"

namespace qtr {

    template<typename T, uint64_t B>
    class MemoryMapVector {
    public:
        static const size_t elementsInBlock = B / sizeof(T);

        static std::unique_ptr<MemoryMapVector<T, B>> create(std::filesystem::path vectorDir);

        static std::unique_ptr<MemoryMapVector<T, B>> open(std::filesystem::path vectorDir);

        const T &operator[](size_t i) const;

        T &operator[](size_t i);

        size_t size() const;

        void push_back(const T &val);

        void push_back(T &&val);

    private:
        MemoryMapVector() = default;

        std::unique_ptr<MemoryMapBytesStorage> _byteStorage;
        size_t _size = 0;
    };

    template<typename T, uint64_t B>
    void MemoryMapVector<T, B>::push_back(T &&val) {
        void *ptr = _byteStorage->allocate(sizeof(T));
        new(ptr) T(std::move(val));
        _size++;
    }

    template<typename T, uint64_t B>
    void MemoryMapVector<T, B>::push_back(const T &val) {
        void *ptr = _byteStorage->allocate(sizeof(T));
        new(ptr) T(val);
        _size++;
    }

    template<typename T, uint64_t B>
    size_t MemoryMapVector<T, B>::size() const {
        return _size;
    }

    template<typename T, uint64_t B>
    T &MemoryMapVector<T, B>::operator[](size_t i) {
        return *const_cast<T *>(const_cast<MemoryMapVector<T, B>*>(this)->operator[](i));

    }

    template<typename T, uint64_t B>
    const T &MemoryMapVector<T, B>::operator[](size_t i) const {
        return *static_cast<const T *>(_byteStorage->getBytes(i / elementsInBlock, i % elementsInBlock));
    }

    template<typename T, uint64_t B>
    std::unique_ptr<MemoryMapVector<T, B>> MemoryMapVector<T, B>::open(std::filesystem::path vectorDir) {
        MemoryMapVector<T, B> *vector = new MemoryMapVector<T, B>;
        vector->_byteStorage = MemoryMapBytesStorage::open(vectorDir);
        vector->_size =
                vector->_byteStorage->filesCount() * elementsInBlock + vector->_byteStorage->lastFileSize() / sizeof(T);
        return std::unique_ptr<MemoryMapVector<T, B>>(vector);
    }

    template<typename T, uint64_t B>
    std::unique_ptr<MemoryMapVector<T, B>>
    MemoryMapVector<T, B>::create(std::filesystem::path vectorDir) {
        MemoryMapVector<T, B> *vector = new MemoryMapVector<T, B>;
        vector->_byteStorage = MemoryMapBytesStorage::create(vectorDir, B);
        vector->_size = 0;
        return std::unique_ptr<MemoryMapVector<T, B>>(vector);
    }

} // qtr
