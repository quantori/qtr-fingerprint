#pragma once

#include <cstdio>
#include <cstring>

namespace qtr {

    template<size_t buf_size>
    class BufferedWriter {
    private:
        FILE *_file;
        char _buf[buf_size];
        size_t _buf_ptr;

    public:
        BufferedWriter(const char *filePath) : _buf_ptr(0) {
            _file = std::fopen(filePath, "w");
        }

        BufferedWriter(const std::filesystem::path &filePath) : BufferedWriter(filePath.c_str()) {}

        BufferedWriter(const std::string &filePath) : BufferedWriter(filePath.c_str()) {}

        void flush() {
            fwrite(_buf, 1, _buf_ptr, _file);
            _buf_ptr = 0;
        }

        BufferedWriter &put(char symbol) {
            if (_buf_ptr == buf_size)
                flush();
            _buf[_buf_ptr++] = symbol;
            return *this;
        }

        BufferedWriter &write(const char *s, size_t count) {
            while (count != 0) {
                if (count <= buf_size - _buf_ptr) {
                    memcpy(_buf + _buf_ptr, s, count);
                    s += count;
                    _buf_ptr += count;
                    count = 0;
                } else {
                    memcpy(_buf + _buf_ptr, s, buf_size - _buf_ptr);
                    s += buf_size - _buf_ptr;
                    count -= buf_size - _buf_ptr;
                    _buf_ptr = buf_size;
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