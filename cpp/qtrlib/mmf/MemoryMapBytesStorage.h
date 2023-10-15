#pragma once

#include <vector>

#include "MemoryMapFile.h"

namespace qtr {

    class MemoryMapBytesStorage {
    public:
        static std::unique_ptr<MemoryMapBytesStorage> create(std::filesystem::path storageDir, uint64_t blockSize);

        static std::unique_ptr<MemoryMapBytesStorage> open(std::filesystem::path storageDir);

        void* allocate(size_t bytesCount);

        void* getBytes(size_t fileId, size_t offset);

        const void* getBytes(size_t fileId, size_t offset) const;

        size_t filesCount() const;

        size_t lastFileSize() const;
    private:
        MemoryMapBytesStorage() = default;

        std::filesystem::path _storageDir;
        std::vector<std::unique_ptr<MemoryMapFile>> _memoryMapFiles;
        uint64_t _filesCount = 0;
        uint64_t _blockSize = 0;
    };

} // qtr
