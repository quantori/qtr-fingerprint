#pragma once

#include <cstdio>
#include <cstring>

namespace qtr {

    static const size_t DefaultReaderBufferSize = 1ull << 13;

    template<size_t bufSize = DefaultReaderBufferSize>
    class BufferedReader {
    private:
        FILE *_file;
        char _buf[bufSize];
        size_t _currentBufSize;
        size_t _bufPtr;
        bool _eof;

    private:
        bool checkEof() {
            if (_bufPtr == _currentBufSize) {
                _bufPtr = 0;
                _currentBufSize = fread(_buf, 1, bufSize, _file);
            }
            return _eof = _currentBufSize == 0;
        }

    public:
        BufferedReader(const char *filePath) : _currentBufSize(0), _bufPtr(0), _eof(false) {
            _file = std::fopen(filePath, "r");
            if (_file == nullptr) {
                LOG(INFO) << "Can't open file " << filePath << ", errno " << errno;
            }
        }

        BufferedReader(const std::filesystem::path &filePath) : BufferedReader(filePath.c_str()) {}

        BufferedReader(const std::string &filePath) : BufferedReader(filePath.c_str()) {}

        bool eof() const {
            return _eof;
        }

        int get() {
            return checkEof() ? EOF : (int) (unsigned char) _buf[_bufPtr++];
        }

        int peek() {
            return checkEof() ? EOF : (int) (unsigned char) _buf[_bufPtr];
        }

        BufferedReader &read(char *s, size_t count) {
            while (count != 0 && !checkEof()) {
                if (count <= _currentBufSize - _bufPtr) {
                    memcpy(s, _buf + _bufPtr, count);
                    _bufPtr += count;
                    s += count;
                    count = 0;
                } else {
                    memcpy(s, _buf + _bufPtr, _currentBufSize - _bufPtr);
                    count -= _currentBufSize - _bufPtr;
                    s += _currentBufSize - _bufPtr;
                    _bufPtr = _currentBufSize;
                }
            }
            return *this;
        }

        ~BufferedReader() {
            fclose(_file);
        }
    };
}