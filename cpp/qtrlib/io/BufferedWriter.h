#pragma once

#include <cstdio>
#include <cstring>

namespace qtr {

    // TODO^ test this class
    template<size_t buf_size>
    class BufferedWriter {
    private:
        FILE *_file;
        char _buf[buf_size];
        size_t _buf_ptr;

    public:
        BufferedWriter(const char *fileName) : _buf_ptr(0) {
            _file = std::fopen(fileName, "w");
        }

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
                    count = 0;
                } else {
                    memcpy(_buf + _buf_ptr, s, buf_size - _buf_ptr);
                    s += buf_size - _buf_ptr;
                    count -= buf_size - _buf_ptr;
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