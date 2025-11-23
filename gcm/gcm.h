#ifndef GCM_H
#define GCM_H


#include <wmmintrin.h>
#include <emmintrin.h>
#include <smmintrin.h>
#include <QString>

/**
 * Reflect a complete 128bit-string. bit0<>bit127, bit1<>bit126, ...
 *
 * First, the string is split into two 128bit temp registers containing the lower (tmp1) and higher (tmp2) nimbles (4bit) of each byte.
 * (Achieved through srli_epi16 and AND-masks)
 *
 * Next, pre-defined HIGHER and LOWER masks are shuffled, based on the values in tmp1 and tmp2.
 * The values of lower nimbles in tmp1 are used to index the values in the HIGHER mask (e.g. for index 0: result[0] = HIGHER-MASK[tmp1[0]]
 * Vice versa for the higher nimbles indexing the values of the LOWER mask.
 * This step uses the _mm_shuffle_epi8 operation and realizes the reflection within each byte.
 *
 * Finally, to reflect the byte order at last, (after xor-ing tmp1 and tmp2), the sum is shuffled according to the BSWAP mask (byte 0<>byte15, byte1<>byte14, ...)
 * @param X the 128bit string
 * @return the reflected string
 */
__m128i reflect_xmm(__m128i X);

void encrypt();
__m128i gfmul(__m128i a, __m128i b);
__m128i gfmul(__m128i a, __m128i b, __m128i q);
__m128i bitshift_left(__m128i a, unsigned char count);
struct int256 bitshift_left256(__m128i i32, __m128i i10, unsigned char count);
__m128i bitshift_right(__m128i a, unsigned char count);
QString print128_hex_lanes(__m128i var);
void gfmul_test();
struct int256{
    __m128i t10;
    __m128i t32;
};
static __m128i Q = _mm_set_epi32(0, 0, 0, 0x00000087);
static __m128i Q_r = _mm_set_epi32(0, 0, 0xc2000000, 0);
static __m128i ZERO = _mm_setzero_si128();

#endif // GCM_H
