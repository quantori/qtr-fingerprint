#pragma once

#include <cstdio>
#include <cstring>

namespace qtr {

    template<size_t buf_size>
    class BufferedReader {
    private:
        FILE *_file;
        char _buf[buf_size];
        size_t _current_buf_size;
        size_t _buf_ptr;
        bool _eof;

    private:
        bool checkEof() {
            if (_buf_ptr == _current_buf_size) {
                _buf_ptr = 0;
                _current_buf_size = fread(_buf, 1, buf_size, _file);
            }
            return _eof = _current_buf_size == 0;
        }

    public:
        BufferedReader(const char *filePath) : _current_buf_size(0), _buf_ptr(0), _eof(false) {
            _file = std::fopen(filePath, "r");
        }

        BufferedReader(const std::filesystem::path &filePath) : BufferedReader(filePath.c_str()) {}

        BufferedReader(const std::string &filePath) : BufferedReader(filePath.c_str()) {}

        bool eof() const {
            return _eof;
        }

        int get() {
            return checkEof() ? EOF : (int) (unsigned char) _buf[_buf_ptr++];
        }

        int peek() {
            return checkEof() ? EOF : (int) (unsigned char) _buf[_buf_ptr];
        }

        BufferedReader &read(char *s, size_t count) {
            while (count != 0 && !checkEof()) {
                if (count <= _current_buf_size - _buf_ptr) {
                    memcpy(s, _buf + _buf_ptr, count);
                    _buf_ptr += count;
                    s += count;
                    count = 0;
                } else {
                    memcpy(s, _buf + _buf_ptr, _current_buf_size - _buf_ptr);
                    count -= _current_buf_size - _buf_ptr;
                    s += _current_buf_size - _buf_ptr;
                    _buf_ptr = _current_buf_size;
                }
            }
            return *this;
        }

        ~BufferedReader() {
            fclose(_file);
        }
    };
}