#pragma once

#include <cstdio>
#include <cstring>

namespace qtr {

    static const size_t DefaultWriterBufferSize = 1ull << 13;

    template<size_t buf_size = DefaultWriterBufferSize>
    class BufferedWriter {
    private:
        FILE *_file;
        char _buf[buf_size];
        size_t _bufPtr;

    public:
        BufferedWriter(const char *filePath) : _bufPtr(0) {
            std::ifstream kek(filePath);
            _file = std::fopen(filePath, "w");
            if (_file == nullptr) {
                LOG(INFO) << "Can't open file " << filePath << ", errno " << errno;
            }
        }

        BufferedWriter(const std::filesystem::path &filePath) : BufferedWriter(filePath.c_str()) {}

        BufferedWriter(const std::string &filePath) : BufferedWriter(filePath.c_str()) {}

        void flush() {
            fwrite(_buf, 1, _bufPtr, _file);
            _bufPtr = 0;
        }

        BufferedWriter &put(char symbol) {
            if (_bufPtr == buf_size)
                flush();
            _buf[_bufPtr++] = symbol;
            return *this;
        }

        BufferedWriter &write(const char *s, size_t count) {
            while (count != 0) {
                if (count <= buf_size - _bufPtr) {
                    memcpy(_buf + _bufPtr, s, count);
                    s += count;
                    _bufPtr += count;
                    count = 0;
                } else {
                    memcpy(_buf + _bufPtr, s, buf_size - _bufPtr);
                    s += buf_size - _bufPtr;
                    count -= buf_size - _bufPtr;
                    _bufPtr = buf_size;
                    flush();
                }
            }
            return *this;
        }

        ~BufferedWriter() {
            flush();
            fclose(_file);
        }
    };

}