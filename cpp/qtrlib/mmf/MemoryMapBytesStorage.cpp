#include "MemoryMapBytesStorage.h"

#include <fstream>

using namespace std;
using namespace qtr;

namespace {
    size_t getFileSize(const filesystem::path &filePath) {
        ifstream file(filePath, ios::binary | ios::ate);
        return file.tellg();
    }

    size_t getFilesCount(const filesystem::path &dirPath) {
        return std::distance(filesystem::directory_iterator(dirPath), filesystem::directory_iterator());
    }
}

std::unique_ptr<MemoryMapBytesStorage>
MemoryMapBytesStorage::create(std::filesystem::path storageDir, uint64_t blockSize) {
    MemoryMapBytesStorage *storage = new MemoryMapBytesStorage();
    storage->_storageDir = std::move(storageDir);
    if (filesystem::exists(storage->_storageDir)) {
        throw std::runtime_error("MemoryMapBytesStorage: storage directory already exists");
    }
    filesystem::create_directory(storage->_storageDir);
    auto firstFile = MemoryMapFile::create(storage->_storageDir / "0", blockSize);
    storage->_memoryMapFiles.push_back(std::move(firstFile));
    storage->_filesCount = 1;
    storage->_blockSize = blockSize;
    return unique_ptr<MemoryMapBytesStorage>(storage);
}

std::unique_ptr<MemoryMapBytesStorage> MemoryMapBytesStorage::open(std::filesystem::path storageDir) {
    filesystem::path firstFilePath = storageDir / "0";
    if (!filesystem::is_regular_file(firstFilePath)) {
        throw std::runtime_error("MemoryMapBytesStorage: storage is damaged or not exits");
    }
    MemoryMapBytesStorage *storage = new MemoryMapBytesStorage();
    storage->_storageDir = std::move(storageDir);
    storage->_blockSize = getFileSize(firstFilePath);
    storage->_filesCount = getFilesCount(storage->_storageDir);
    auto firstFile = MemoryMapFile::open(storage->_storageDir / "0");
    storage->_memoryMapFiles.push_back(std::move(firstFile));
    for (uint64_t i = 1; i < storage->_filesCount; i++) {
        storage->_memoryMapFiles.push_back(MemoryMapFile::open(storage->_storageDir / to_string(i)));
    }
    return unique_ptr<MemoryMapBytesStorage>(storage);
}

void *MemoryMapBytesStorage::allocate(size_t bytesCount) {
    void *ptr = _memoryMapFiles.back()->allocate(bytesCount);
    if (ptr != nullptr)
        return ptr;
    _memoryMapFiles.push_back(MemoryMapFile::create(_storageDir / to_string(_memoryMapFiles.size()), _blockSize));
    return _memoryMapFiles.back()->allocate(bytesCount);
}

void *MemoryMapBytesStorage::getBytes(size_t fileId, size_t offset) {
    return _memoryMapFiles[fileId]->ptr(offset);
}

const void *MemoryMapBytesStorage::getBytes(size_t fileId, size_t offset) const {
    return _memoryMapFiles[fileId]->ptr(offset);
}

size_t MemoryMapBytesStorage::filesCount() const {
    return _filesCount;
}

size_t MemoryMapBytesStorage::lastFileSize() const {
    return _memoryMapFiles.back()->size();
}
