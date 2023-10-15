#include "MemoryMapFile.h"

#include <fcntl.h>
#include <cstring>
#include <unistd.h>
#include <sys/mman.h>

using namespace qtr;
using namespace std;

namespace {
    int openMemoryMapFile(const std::filesystem::path &filePath) {
        int fd;
        if ((fd = open(filePath.c_str(), O_RDWR | O_CREAT, S_IRUSR | S_IWUSR)) == -1) {
            throw std::runtime_error("MemoryMapFile: Cannot open file: " + string(strerror(errno)));
        }
        return fd;
    }

    void *memoryMap(size_t len, int fd) {
        void *ptr = mmap(NULL, len, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
        if (ptr == MAP_FAILED) {
            throw std::runtime_error("MemoryMapFile: Cannot map the file: " + string(strerror(errno)));
        }
        return (char *) ptr + sizeof(MemoryMapFile::Header);
    }

    uint64_t fileSize(uint64_t capacity) {
        return sizeof(MemoryMapFile::Header) + capacity;
    }
}

MemoryMapFile::~MemoryMapFile() {
    if (_ptr != nullptr) {
        munmap(_ptr, capacity());
    }

    if (_fd != -1) {
        close(_fd);
    }
}

void *MemoryMapFile::ptr(ptrdiff_t offset) {
    return (char *) _ptr + offset;
}

const void *MemoryMapFile::ptr(ptrdiff_t offset) const {
    return (char *) _ptr + offset;
}


void *MemoryMapFile::allocate(size_t bytesCount) {
    if (size() + bytesCount > _capacity())
        return nullptr;
    _size() += bytesCount;
    return ptr(size() - bytesCount);
}

const MemoryMapFile::Header &MemoryMapFile::header() const {
    return *(static_cast<const Header *>(ptr()) - 1);
}

MemoryMapFile::Header &MemoryMapFile::_header() {
    return *(static_cast<Header *>(ptr()) - 1);
}

std::unique_ptr<MemoryMapFile> MemoryMapFile::create(const std::filesystem::path &filePath, uint64_t capacity) {
    if (filesystem::exists(filePath)) {
        throw std::runtime_error("Cannot create memory map file " + filePath.string() + " because it already exists");
    }
    MemoryMapFile *mmf = new MemoryMapFile();
    uint64_t maxFileSize = fileSize(capacity);
    mmf->_fd = openMemoryMapFile(filePath);
    if (ftruncate(mmf->_fd, maxFileSize) < 0) {
        throw std::runtime_error("MemoryMapFile: Cannot truncate file to size=" + to_string(maxFileSize));
    }
    mmf->_ptr = memoryMap(maxFileSize, mmf->_fd);
    mmf->_capacity() = capacity;
    mmf->_size() = 0;
    return unique_ptr<MemoryMapFile>(mmf);
}

std::unique_ptr<MemoryMapFile> MemoryMapFile::open(const std::filesystem::path &filePath) {
    if (!filesystem::is_regular_file(filePath)) {
        throw std::runtime_error(
                "Cannot open memory map file " + filePath.string() +
                " because it does not exist or not a regular file");
    }
    MemoryMapFile *mmf = new MemoryMapFile();
    mmf->_fd = openMemoryMapFile(filePath);
    Header header;
    read(mmf->_fd, &header, sizeof(header));
    mmf->_ptr = memoryMap(fileSize(header.capacity), mmf->_fd);
    return unique_ptr<MemoryMapFile>(mmf);
}

uint64_t &MemoryMapFile::_capacity() {
    return _header().capacity;
}

uint64_t MemoryMapFile::capacity() const {
    return header().capacity;
}

uint64_t &MemoryMapFile::_size() {
    return _header().size;
}

uint64_t MemoryMapFile::size() const {
    return header().size;
}


