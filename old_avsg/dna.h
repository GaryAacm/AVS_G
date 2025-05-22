#ifndef AVS_API_DNA_H
#define AVS_API_DNA_H

// cross-platform byte swap
#ifdef _MSC_VER

#include <stdlib.h>
#define bswap_32(x) _byteswap_ulong(x)
#define bswap_64(x) _byteswap_uint64(x)

#elif defined(__APPLE__)

// Mac OS X / Darwin features
#include <libkern/OSByteOrder.h>
#define bswap_32(x) OSSwapInt32(x)
#define bswap_64(x) OSSwapInt64(x)

#elif defined(__sun) || defined(sun)

#include <sys/byteorder.h>
#define bswap_32(x) BSWAP_32(x)
#define bswap_64(x) BSWAP_64(x)

#elif defined(__FreeBSD__)

#include <sys/endian.h>
#define bswap_32(x) bswap32(x)
#define bswap_64(x) bswap64(x)

#elif defined(__OpenBSD__)

#include <sys/types.h>
#define bswap_32(x) swap32(x)
#define bswap_64(x) swap64(x)

#elif defined(__NetBSD__)

#include <sys/types.h>
#include <machine/bswap.h>
#if defined(__BSWAP_RENAME) && !defined(__bswap_32)
#define bswap_32(x) bswap32(x)
#define bswap_64(x) bswap64(x)
#endif

#else

#include <byteswap.h>

#endif

#include <cstdint>
#include <algorithm>
#include "string-view-lite/include/nonstd/string_view.hpp"
#include <string>
#include <vector>

using u8 = uint8_t;
using u16 = uint16_t;
using u32 = uint32_t;
using u64 = uint64_t;

//constexpr decltype(auto) base_encoder_table() {
//    std::array<u8, 256> table{};
//    for (size_t i = 0; i < 256; i++) {
//        switch (i) {
//            case 'A':
//            case 'a':
//                table[i] = 0; // 00
//                break;
//            case 'C':
//            case 'c':
//                table[i] = 1; // 01
//                break;
//            case 'G':
//            case 'g':
//                table[i] = 2; // 10
//                break;
//            case 'T':
//            case 't':
//                table[i] = 3; // 11
//                break;
//            default: // 其他非碱基的字符一律编码为4(比如序列中经常出现的'N')
//                table[i] = 4;
//                break;
//        }
//    }
//    return table;
//}

static constexpr uint8_t encode_table[256] = { // A : 0; C : 1; G : 2; T : 3
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

//static constexpr const auto encode_table = base_encoder_table();
static constexpr const char base_decoder_table[] = "ACGTN-acgtn*";

inline u8 byte_2bit_unit_reverse_8(u8 number) {
    number = ((number & 0x33ULL) << 2ULL) | ((number & 0xCCULL) >> 2ULL);
    number = ((number & 0x0FULL) << 4ULL) | ((number & 0xF0ULL) >> 4ULL);
    return number;
}

inline u16 byte_2bit_unit_reverse_16(u16 number) {
    number = ((number & 0x3333ULL) << 2ULL) | ((number & 0xCCCCULL) >> 2ULL);
    number = ((number & 0x0F0FULL) << 4ULL) | ((number & 0xF0F0ULL) >> 4ULL);
    return number;
}

inline u32 byte_2bit_unit_reverse_32(u32 number) {
    number = ((number & 0x33333333ULL) << 2ULL) | ((number & 0xCCCCCCCCULL) >> 2ULL);
    number = ((number & 0x0F0F0F0FULL) << 4ULL) | ((number & 0xF0F0F0F0ULL) >> 4ULL);
    return number;
}

inline u64 byte_2bit_unit_reverse_64(u64 number) {
    number = ((number & 0x3333333333333333ULL) << 2ULL) | ((number & 0xCCCCCCCCCCCCCCCCULL) >> 2ULL);
    number = ((number & 0x0F0F0F0F0F0F0F0FULL) << 4ULL) | ((number & 0xF0F0F0F0F0F0F0F0ULL) >> 4ULL);
    return number;
}

inline u8 byte_1bit_unit_reverse_8(u8 number) {
    number = ((number & 0x55ULL) << 1ULL) | ((number & 0xAAULL) >> 1ULL);
    return byte_2bit_unit_reverse_8(number);
}

inline u16 byte_1bit_unit_reverse_16(u16 number) {
    number = ((number & 0x5555ULL) << 1ULL) | ((number & 0xAAAAULL) >> 1ULL);
    return byte_2bit_unit_reverse_16(number);
}

inline u32 byte_1bit_unit_reverse_32(u32 number) {
    number = ((number & 0x55555555ULL) << 1ULL) | ((number & 0xAAAAAAAAULL) >> 1ULL);
    return byte_2bit_unit_reverse_32(number);
}

inline u64 byte_1bit_unit_reverse_64(u64 number) {
    number = ((number & 0x5555555555555555ULL) << 1ULL) | ((number & 0xAAAAAAAAAAAAAAAAULL) >> 1ULL);
    return byte_2bit_unit_reverse_64(number);
}

static inline uint64_t get_reverse_kmer(uint64_t kmer, uint8_t kmer_size) {
    kmer = ~kmer;
    kmer = ((kmer & 0x3333333333333333ULL) << 2ULL) | ((kmer & 0xCCCCCCCCCCCCCCCCULL) >> 2ULL);
    kmer = ((kmer & 0x0F0F0F0F0F0F0F0FULL) << 4ULL) | ((kmer & 0xF0F0F0F0F0F0F0F0ULL) >> 4ULL);
    kmer = bswap_64(kmer);
    return (kmer >> ((32ULL - kmer_size) << 1ULL));
}

static inline uint64_t kmer_to_u64(const nonstd::string_view &kmer) {
    // assert(kmer.size()<=32);
    uint64_t bits = 0, c;
    u8 idx2 = 0;
    size_t kmer_bits = kmer.size() * 2;
    for (uint8_t ch : kmer) {
        c = encode_table[ch] & 0x03U;
        bits |= (c << (kmer_bits - 2 - idx2));
        idx2 = (idx2 + 2U) & 0x3FU;
    }
    return bits;
}

#endif