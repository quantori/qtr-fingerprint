#pragma once

#include <vector>
#include "base_cpp/array.h"

template<typename T>
class IndigoQueryFingerprint {
public:
    explicit IndigoQueryFingerprint(const indigo::Array<T> &fingerprint) {
        for (int i = 0; i < fingerprint.size(); i++) {
            byte b = fingerprint.at(i);
            if (b == 0) {
                continue;
            }
            _items.emplace_back(i, b);
        }
    }

    [[nodiscard]] bool isSubFingerprint(const indigo::Array<T> &fingerprint) const {
        const T* fp = fingerprint.ptr();
        for (auto &[idx, b]: _items) {
            byte val = fp[idx];
            if ((b & val) != b) {
                return false;
            }
        }
        return true;
    }

private:
    std::vector<std::pair<int, T>> _items;
};

//template<typename T>
//class IndigoQueryFingerprint {
//public:
//    inline explicit IndigoQueryFingerprint(const indigo::Array<byte> &fingerprint) {
//        assert(fingerprint.size() % sizeof(uint32_t) == 0);
//        _fp.copy(reinterpret_cast<const uint32_t *>(fingerprint.ptr()), fingerprint.size() / sizeof(uint32_t));
//    }
//
//    [[nodiscard]] inline bool isSubFingerprint(const indigo::Array<byte> &fingerprint) const {
//        const uint32_t *fp1 = _fp.ptr();
//        const uint32_t *fp2 = reinterpret_cast<const uint32_t *>(fingerprint.ptr());
//        const size_t length = _fp.size();
//
//        constexpr std::size_t words_per_vec = 256 / (8 * sizeof(uint32_t));   // 8 for u32
//        std::size_t i = 0;
//
//        for (; i + words_per_vec <= length; i += words_per_vec) {
//            __m256i v1 = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(fp1 + i));
//            __m256i v2 = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(fp2 + i));
//
//            __m256i vand = _mm256_and_si256(v1, v2); // (fp1 & fp2)
//            __m256i cmp = _mm256_cmpeq_epi32(vand, v1); // compare 8×u32
//
//            if (_mm256_movemask_epi8(cmp) != -1)
//                return false;
//        }
//
//        for (; i < length; ++i)
//            if ((fp1[i] & fp2[i]) != fp1[i])
//                return false;
//
//        return true;
//    }
//
//private:
//    indigo::Array<uint32_t> _fp;
//};

