#pragma once

#include <filesystem>

namespace qtr {

    class MemoryMapFile {
    public:
        struct Header {
            uint64_t capacity;
            uint64_t size;
        };

        static MemoryMapFile create(const std::filesystem::path &filePath, uint64_t capacity);

        static MemoryMapFile open(const std::filesystem::path &filePath);

        ~MemoryMapFile();

        MemoryMapFile &operator=(const MemoryMapFile &) = delete;

        MemoryMapFile(const MemoryMapFile &) = delete;

        MemoryMapFile &operator=(MemoryMapFile &&) = delete;

        MemoryMapFile(MemoryMapFile &&) = default;

        void *ptr(ptrdiff_t offset = 0);

        const void *ptr(ptrdiff_t offset = 0) const;

        void *allocate(size_t bytesCount);

        const Header &header() const;

        uint64_t capacity() const;

        uint64_t size() const;

    private:
        Header &_header();

        uint64_t &_capacity();

        uint64_t &_size();

        MemoryMapFile() = default;

        void *_ptr = nullptr;
        int _fd = -1;
    };
} // qtr

