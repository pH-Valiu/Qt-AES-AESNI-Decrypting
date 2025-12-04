#include <gcm/gfmul.h>
#include "gcm/benchmarkutil.h"
#include <QDebug>

extern "C" void gfmul_reflected_avx512(__m512i* a, __m512i* b, __m512i* result);
extern "C" void gfmul_reflected_avx128(__m128i* a, __m128i* b, __m128i* result);

/**
 * @brief reflect_xmm code from Intel Doc. A
 * @param X
 * @return
 */
__m128i reflect_xmm(__m128i X){
    __m128i tmp1,tmp2;
    __m128i AND_MASK =
        _mm_set_epi32(0x0f0f0f0f, 0x0f0f0f0f, 0x0f0f0f0f, 0x0f0f0f0f);
    __m128i LOWER_MASK =
        _mm_set_epi32(0x0f070b03, 0x0d050901, 0x0e060a02, 0x0c040800);
    __m128i HIGHER_MASK =
        _mm_set_epi32(0xf070b030, 0xd0509010, 0xe060a020, 0xc0408000);
    __m128i BSWAP_MASK =
        _mm_set_epi8(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15);
    tmp2 = _mm_srli_epi16(X, 4);        // splits the string into 8x 16bit words (e.g. [0xABCD]). Each is >> 4, resulting in [0x0ABC]. This masked with [0x0f0f], gives [0x0A0C] directly resembling the higher nimbles.
    tmp1 = _mm_and_si128(X, AND_MASK);
    tmp2 = _mm_and_si128(tmp2, AND_MASK);
    tmp1 = _mm_shuffle_epi8(HIGHER_MASK ,tmp1);
    tmp2 = _mm_shuffle_epi8(LOWER_MASK ,tmp2);
    tmp1 = _mm_xor_si128(tmp1, tmp2);
    return _mm_shuffle_epi8(tmp1, BSWAP_MASK);
};

__m128i gfmul(const __m128i& a, const __m128i& b){
    // Step 1: Multiply
    __m128i a0b0 = _mm_clmulepi64_si128(a, b, 0x00);
    __m128i a0b1 = _mm_clmulepi64_si128(a, b, 0x10);
    __m128i a1b0 = _mm_clmulepi64_si128(a, b, 0x01);
    __m128i a1b1 = _mm_clmulepi64_si128(a, b, 0x11);

    __m128i mid = _mm_xor_si128(a0b1, a1b0); // computes mid = A0B1 + A1B0

    __m128i c01 = _mm_xor_si128(a0b0, _mm_slli_si128(mid, 8)); // computes C[1:0] = A0B0 + (mid << x^64)
    __m128i c23 = _mm_xor_si128(a1b1, _mm_srli_si128(mid, 8)); // computes C[3:2] = A1B1 + (mid >> x^64)

    // Step 2.1: Reduce
    __m128i x = _mm_clmulepi64_si128(a1b1, Q, 0x01);  // computes C[3] * Q = upper(A1B1) * lower(Q)
    c01 = _mm_xor_si128(c01, _mm_slli_si128(x, 8)); // add lower half of x (X[0]) to upper part of C[1:0]
    c23 = _mm_xor_si128(c23, _mm_srli_si128(x, 8)); // add higher half of x (X[1]) to lower part of C[3:2] (higher part is just dangling around)

    // Step 2.2: Reduce
    x = _mm_clmulepi64_si128(c23, Q, 0x00);         // works because higher part is not used in the calculation
    c01 = _mm_xor_si128(c01, x);

    return c01;
}


/*
 * Please note: We use the Symbol "GFMUL" to refer to "Carry-Less-Multiply-And-Reduce",
 *  while Intel Doc. B uses "*" to refer to "Carry-Less-Multiply-And-Reduce".
 * Meanwhile, we use "CLMUL" to refer to just "Carry-Less-Multiply",
 *  while Intel Doc. B uses "GFMUL64" to refer to just "Carry-Less-Multiply".
 * --> Note the similarity but difference between "GFMUL" and "GFMUL64".
 *
 * Meanwhile, Intel Doc. A mostly uses "\cdot" to refer to "Carry-Less-Multiply" and ("*" or "GFMUL128") to refer to "Carry-Less-Multiply-And-Reduce"
 */
/**
 * @brief gfmul_k_optimized This should only be used when we have b input directly from powers of H.
 * So, as a standalone, internally, another gfmul() call is issued. --> Per se: Bad
 * @param a
 * @param b
 * @return
 */
__m128i gfmul_k_optimized(const __m128i& a, const __m128i& b){
    // Step 0: Pre-requesites
    //__m128i k = _mm_xor_si128(_mm_clmulepi64_si128(b, q, 0x01), _mm_slli_si128(b, 8)); This line achieved the same as the one below (in reversed case think about the bitshift which could result in greater >128 bit intermediate result)
    __m128i k = _mm_xor_si128(gfmul(_mm_srli_si128(b, 8), Q), _mm_slli_si128(b, 8));     // K = GFMUL(B[1], Q) + B[0]*x^64

    // Step 1: Multiply
    __m128i a0b0 = _mm_clmulepi64_si128(a, b, 0x00);
    __m128i a0b1 = _mm_clmulepi64_si128(a, b, 0x10);
    __m128i a1k0 = _mm_clmulepi64_si128(a, k, 0x01);
    __m128i a1k1 = _mm_clmulepi64_si128(a, k, 0x11);

    __m128i lower  = _mm_xor_si128(a0b0, a1k0);     // computes lower  = a0b0 + a1k0
    __m128i higher = _mm_xor_si128(a0b1, a1k1);     // computes higher = a0b1 + a1k1

    __m128i c01 = _mm_xor_si128(lower, _mm_slli_si128(higher, 8));  // c01 = lower + (higher << 64)
    __m128i c23 = _mm_srli_si128(higher, 8);                        // c02 = higher >> 64           (highest 64 bits are 0)

    // Step 2: Reduce
    __m128i x = _mm_clmulepi64_si128(c23, Q, 0x00);
    c01 = _mm_xor_si128(x, c01);

    return c01;
}



/**
 * @brief gfmul_reflected
 * This approach mostly strictly follows the "Simple Bit-Reflected GHASH Algorithm" of Table 9 in Doc B.
 * Important note: while the table says that inputs a and b are reflected, they are indeed not.
 * Additionally, after the CLMUL(A,B), the result must be shifted by one to the left (CLMUL(A,B) << 1) to accomadate for the fact that we are intrinsically working on reflected values.
 * Meanwhile, the reduction polynom is correctly Q' = reflect_64(Q>>1) = epi32(0...0, 0...0, 0xc2000000, 0...0)
 * The returnted value must not be reflected but instead immediately represents the correct output as if we had done:
 * OUT = (GFMUL(A', B'))' = (OUT')' = OUT
 *
 * Bear in mind that (CLMUL(A<<1,B))' != (CLMUL(A,B)<<1)'
 *
 * Steps are:
 *
 * 1. Compute CLMUL(A, B) << 1
 * 2. Reduce with Q' (treated as 64bit) following the Doc. B approach
 * @param a
 * @param b
 * @return
 */
__m128i gfmul_reflected(const __m128i& a, const __m128i& b){
    //__m128i Q_r = _mm_set_epi32(0, 0, 0xc2000000, 0);

    // Step 1: Multiply
    __m128i a0b0 = _mm_clmulepi64_si128(a, b, 0x00);
    __m128i a0b1 = _mm_clmulepi64_si128(a, b, 0x10);
    __m128i a1b0 = _mm_clmulepi64_si128(a, b, 0x01);
    __m128i a1b1 = _mm_clmulepi64_si128(a, b, 0x11);

    __m128i mid = _mm_xor_si128(a0b1, a1b0);      // computes mid = A0B1 + A1B0

    __m128i c01 = _mm_xor_si128(a0b0, _mm_slli_si128(mid, 8));    // computes C[1:0] = A0B0 + (mid << x^64)
    __m128i c23 = _mm_xor_si128(a1b1, _mm_srli_si128(mid, 8));    // computes C[3:2] = A1B1 + (mid >> x^64)

    // Step 1.1: Bitshift << 1
    __m128i C01_32bitMSB, C23_32bitMSB, C01_128bitMSB;
    C01_32bitMSB = _mm_srli_epi32(c01, 31);                 // isolate MSB of each 32bit word in C[1:0] and put at index 0 (-31) in each 32bit word
    C23_32bitMSB = _mm_srli_epi32(c23, 31);                 // isolate MSB of each 32bit word in C[3:2] and put at index 0 (-31) in each 32bit word
    c01 = _mm_slli_epi32(c01, 1);                           // shift each 32bit word in C[1:0] by 1 bit to the left, making space at index 0 in each 32bit word
    c23 = _mm_slli_epi32(c23, 1);                           // shift each 32bit word in C[3:2] by 1 bit to the left, making space at index 0 in each 32bit word
    C01_128bitMSB = _mm_srli_si128(C01_32bitMSB, 12);       // bring the all-highest MSB of C[1:0] at index 0 (-127) (this is the carry-over of C[1:0])
    C01_32bitMSB = _mm_slli_si128(C01_32bitMSB, 4);         // rotate the prior isolated MSBs of each 32bit word in C[1:0] by one 32bit block to the left. (the lowest 32bit word is now fully cleared)
    C23_32bitMSB = _mm_slli_si128(C23_32bitMSB, 4);         // rotate the prior isolated MSBs of each 32bit word in C[3:2] by one 32bit block to the left. (the lowest 32bit word is now fully cleared)
    c01 = _mm_or_si128(c01, C01_32bitMSB);                  // add the prior rotated carry-over of each 32bit word onto the next 32bit word (this applies the 32bit-word carry-overs inside C[1:0])
    c23 = _mm_or_si128(c23, C23_32bitMSB);                  // add the prior rotated carry-over of each 32bit word onto the next 32bit word (this applies the 32bit-word carry-overs inside C[1:0])
    c23 = _mm_or_si128(c23, C01_128bitMSB);                 // add the singled out all-highest MSB of C[1:0] into C[3:2] thereby applying the 128bit carry-over from C[1:0] into C[3:2]



    // Step 2.1: Reduce
    __m128i x = _mm_clmulepi64_si128(c01, Q_r, 0x00);
    c23 = _mm_xor_si128(c23, _mm_srli_si128(x, 8));       // add higher half of x (X[1]) to lower part of C[3:2]
    c23 = _mm_xor_si128(c23, _mm_unpacklo_epi64(c01, ZERO));    // add only C[0] to lower part of C[3:2]  (zeroing out C[1])

    c01 = _mm_xor_si128(c01, _mm_slli_si128(x, 8));       // add lower half of x (X[0]) to upper part of C[1:0]

    // Step 2.2: Reduce
    x = _mm_clmulepi64_si128(c01, Q_r, 0x01);               // computes C[1] * Q = higher(C[1:0]) * lower(Q)
    c23 = _mm_xor_si128(c23, x);                          // add full X on C[3:2]
    c23 = _mm_xor_si128(c23, _mm_unpackhi_epi64(ZERO, c01)); // add C[1] on higher C[3:2] (zeroing out C[0])

    return c23;
}

/**
 * It follows the Karatsuba multiplication approach,
 * where we start with [A1:A0] * [B1:B0], and then compute
 * [C1:C0] = A1 * B1
 * [D1:D0] = A0 * B0
 * [E1:E0] = (A1 + A0) * (B1 + B0)
 * Res[3:0] = [C1 : C0 + C1 + E1 + D1 : C0 + E0 + D0 + D1 : D0]
 *                        |    |   |     |    |    |
 * helper sum [F1:F0] = [C1 + D1 + E1 : C0 + D0 + E0]
 *
 * After having computed the full 255 bit result, we apply the same routine as in gfmul_reflected
 * Say, out 255bit result is [G3:G2:G1:G0] (in code this is [C3:C2:C1:C0] just out of a naming thing, is not the same as [C1:C0])
 * then:
 * 1. [G3:G2:G1:G0] << 1
 * 2. Reduce with [0:0:c2000000:0]...
 *
 * Say GHASH is defined over variables Zi as inputs, and Xi as outputs with H as the GHASH key
 * Then the formula for X4 is: (X0 is always 0)
 * X4 =  (X3 + Z4) * H
 *    = ((X2 + Z3) * H) + Z4) * H
 *    = (((((((X0 + Z1) * H) + Z2) * H) + Z3) * H) + Z4) * H
 *    = ((X0 + Z1) * H^4) + (Z2 * H^3) + (Z3 * H^2) + (Z4 * H^1)
 *
 * Now due to the internal nomenclation of the algorithm internally:
 * Y1 is multiplied with H^1
 * Y2 is multiplied with H^2
 * Y3 is multiplied with H^3
 * Y4 is multiplied with H^4
 *
 * Therefore, when invoking this method, the formula:
 * - Z4 corresponds to parameter Y1
 * - Z3 corresponds to parameter Y2
 * - Z2 corresponds to parameter Y3
 * - Z1 corresponds to parameter Y4
 *
 * Hence a typical invocation might look like:
 * ```
 * Z1 = _mm_xor_si128(X, Z1);
 * __m128i X = gfmul_times_four_reflected(Z4, Z3, Z2, Z1, H, H2, H3, H4);
 * ```
 *
 * @brief gfmul_times_four_reflected
 * @param Y1 gets multiplied with H1 (Y1 := Z4)
 * @param Y2 gets multiplied with H2 (Y2 := Z3)
 * @param Y3 gets multiplied with H3 (Y3 := Z2)
 * @param Y4 gets multiplied with H4 (Y4 := Z1)
 * @param H1 H^1
 * @param H2 H^2
 * @param H3 H^3
 * @param H4 H^4
 * @return
 */
__m128i gfmul_times_four_reflected(
    const __m128i &Y1, const __m128i &Y2, const __m128i &Y3, const __m128i &Y4,
    const __m128i &H1, const __m128i &H2, const __m128i &H3, const __m128i &H4)
{
    /*algorithm by Krzysztof Jankowski, Pierre Laurent - Intel*/
    __m128i H1_Y1_lo, H1_Y1_hi,
        H2_Y2_lo, H2_Y2_hi,
        H3_Y3_lo, H3_Y3_hi,
        H4_Y4_lo, H4_Y4_hi,
        lo, hi;
    __m128i tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, c01, c23;

    H1_Y1_lo = _mm_clmulepi64_si128(H1, Y1, 0x00);
    H2_Y2_lo = _mm_clmulepi64_si128(H2, Y2, 0x00);
    H3_Y3_lo = _mm_clmulepi64_si128(H3, Y3, 0x00);
    H4_Y4_lo = _mm_clmulepi64_si128(H4, Y4, 0x00);
    lo = _mm_xor_si128(H1_Y1_lo, H2_Y2_lo);
    lo = _mm_xor_si128(lo, H3_Y3_lo);
    lo = _mm_xor_si128(lo, H4_Y4_lo);
    // --> lo now contains H1Y1_lo + H2Y2_lo + H3Y3_lo + H4Y4_lo

    H1_Y1_hi = _mm_clmulepi64_si128(H1, Y1, 0x11);
    H2_Y2_hi = _mm_clmulepi64_si128(H2, Y2, 0x11);
    H3_Y3_hi = _mm_clmulepi64_si128(H3, Y3, 0x11);
    H4_Y4_hi = _mm_clmulepi64_si128(H4, Y4, 0x11);
    hi = _mm_xor_si128(H1_Y1_hi, H2_Y2_hi);
    hi = _mm_xor_si128(hi, H3_Y3_hi);
    hi = _mm_xor_si128(hi, H4_Y4_hi);
    // --> hi now contains H1Y1_hi + H2Y2_hi + H3Y3_hi + H4Y4_hi

    // compute Hi+Lo (stored in 64bit) for H1 and Y1    (middle term for karatsuba)
    tmp0 = _mm_shuffle_epi32(H1, 78);   // does 64bit halves swap Hi|Lo -> Lo|Hi
    tmp4 = _mm_shuffle_epi32(Y1, 78);
    tmp0 = _mm_xor_si128(tmp0, H1);     // now, each 64 bit half, contains Hi+Lo | Lo+Hi
    tmp4 = _mm_xor_si128(tmp4, Y1);

    // compute Hi+Lo (stored in 64bit) for H2 and Y2    (middle term for karatsuba)
    tmp1 = _mm_shuffle_epi32(H2, 78);
    tmp5 = _mm_shuffle_epi32(Y2, 78);
    tmp1 = _mm_xor_si128(tmp1, H2);
    tmp5 = _mm_xor_si128(tmp5, Y2);

    // compute Hi+Lo (stored in 64bit) for H3 and Y3    (middle term for karatsuba)
    tmp2 = _mm_shuffle_epi32(H3, 78);
    tmp6 = _mm_shuffle_epi32(Y3, 78);
    tmp2 = _mm_xor_si128(tmp2, H3);
    tmp6 = _mm_xor_si128(tmp6, Y3);

    // compute Hi+Lo (stored in 64bit) for H4 and Y4    (middle term for karatsuba)
    tmp3 = _mm_shuffle_epi32(H4, 78);
    tmp7 = _mm_shuffle_epi32(Y4, 78);
    tmp3 = _mm_xor_si128(tmp3, H4);
    tmp7 = _mm_xor_si128(tmp7, Y4);

    // this computes the full middle term for each H0<>Y0 pair (using just the low64bit block on each is fine, we could also use 0x10, 0x01, 0x11)
    tmp0 = _mm_clmulepi64_si128(tmp0, tmp4, 0x00);
    tmp1 = _mm_clmulepi64_si128(tmp1, tmp5, 0x00);
    tmp2 = _mm_clmulepi64_si128(tmp2, tmp6, 0x00);
    tmp3 = _mm_clmulepi64_si128(tmp3, tmp7, 0x00);

    // this computes [F1:F0] = [C1:C0] + [D1:D0] + [E1:E0] (combined over all)
    tmp0 = _mm_xor_si128(tmp0, lo);
    tmp0 = _mm_xor_si128(tmp0, hi);
    tmp0 = _mm_xor_si128(tmp1, tmp0);
    tmp0 = _mm_xor_si128(tmp2, tmp0);
    tmp0 = _mm_xor_si128(tmp3, tmp0);

    tmp4 = _mm_slli_si128(tmp0, 8);     // [F0:00]
    tmp0 = _mm_srli_si128(tmp0, 8);     // [00:F1]
    lo = _mm_xor_si128(tmp4, lo);       // [D1:D0] + [F0:00] = [D1+C0+D0+E0 : D0]
    hi = _mm_xor_si128(tmp0, hi);       // [C1:C0] + [00:F1] = [C1 : C0+C1+D1+E1]

    // now, hi = [C1:C0+C1+D1+E1] || lo = [C0+D0+E0+D1:D0] contains full produt
    c01 = lo;
    c23 = hi;

    // now, we have to reduce again (just copied the reduce operation from gfmul_reflected)
    // therefore, the usage of c01 and c23 does not refer to [C1:C0] but rather refers to
    // the variables [hi:lo] = [C3:C2:C1:C0]

    // Step 1.1: Bitshift << 1
    __m128i C01_32bitMSB, C23_32bitMSB, C01_128bitMSB;
    C01_32bitMSB = _mm_srli_epi32(c01, 31);                 // isolate MSB of each 32bit word in C[1:0] and put at index 0 (-31) in each 32bit word
    C23_32bitMSB = _mm_srli_epi32(c23, 31);                 // isolate MSB of each 32bit word in C[3:2] and put at index 0 (-31) in each 32bit word
    c01 = _mm_slli_epi32(c01, 1);                           // shift each 32bit word in C[1:0] by 1 bit to the left, making space at index 0 in each 32bit word
    c23 = _mm_slli_epi32(c23, 1);                           // shift each 32bit word in C[3:2] by 1 bit to the left, making space at index 0 in each 32bit word
    C01_128bitMSB = _mm_srli_si128(C01_32bitMSB, 12);       // bring the all-highest MSB of C[1:0] at index 0 (-127) (this is the carry-over of C[1:0])
    C01_32bitMSB = _mm_slli_si128(C01_32bitMSB, 4);         // rotate the prior isolated MSBs of each 32bit word in C[1:0] by one 32bit block to the left. (the lowest 32bit word is now fully cleared)
    C23_32bitMSB = _mm_slli_si128(C23_32bitMSB, 4);         // rotate the prior isolated MSBs of each 32bit word in C[3:2] by one 32bit block to the left. (the lowest 32bit word is now fully cleared)
    c01 = _mm_or_si128(c01, C01_32bitMSB);                  // add the prior rotated carry-over of each 32bit word onto the next 32bit word (this applies the 32bit-word carry-overs inside C[1:0])
    c23 = _mm_or_si128(c23, C23_32bitMSB);                  // add the prior rotated carry-over of each 32bit word onto the next 32bit word (this applies the 32bit-word carry-overs inside C[1:0])
    c23 = _mm_or_si128(c23, C01_128bitMSB);                 // add the singled out all-highest MSB of C[1:0] into C[3:2] thereby applying the 128bit carry-over from C[1:0] into C[3:2]



    // Step 2.1: Reduce
    __m128i x = _mm_clmulepi64_si128(c01, Q_r, 0x00);       // computes C[0] * Q = lower(C[1:0]) * lower(Q)
    c23 = _mm_xor_si128(c23, _mm_srli_si128(x, 8));       // add higher half of x (X[1]) to lower part of C[3:2]
    c23 = _mm_xor_si128(c23, _mm_unpacklo_epi64(c01, ZERO));    // add only C[0] to lower part of C[3:2]  (zeroing out C[1])

    c01 = _mm_xor_si128(c01, _mm_slli_si128(x, 8));       // add lower half of x (X[0]) to upper part of C[1:0]

    // Step 2.2: Reduce
    x = _mm_clmulepi64_si128(c01, Q_r, 0x01);               // computes C[1] * Q = higher(C[1:0]) * lower(Q)
    c23 = _mm_xor_si128(c23, x);                          // add full X on C[3:2]
    c23 = _mm_xor_si128(c23, _mm_unpackhi_epi64(ZERO, c01)); // add C[1] on higher C[3:2] (zeroing out C[0])

    return c23;
}


#ifdef AVX512_SUPPORT
/**
 * @brief gfmul_reflected_avx512_parallel
 * This approach mostly strictly follows the "Simple Bit-Reflected GHASH Algorithm" of Table 9 in Doc B.
 * But with the addition of operating on four quadwords (__m128i) at the same time.
 * One call of this method is equivalent to four gfmul_reflected ones where each time a different pair is passed.
 *
 * Important note: while the table says that inputs a and b are reflected, they are indeed not.
 * Additionally, after the CLMUL(A,B), the result must be shifted by one to the left (CLMUL(A,B) << 1) to accomadate for the fact that we are intrinsically working on reflected values.
 * Meanwhile, the reduction polynom is correctly Q' = reflect_64(Q>>1) = epi32(0...0, 0...0, 0xc2000000, 0...0)
 * The returnted value must not be reflected but instead immediately represents the correct output as if we had done:
 * OUT = (GFMUL(A', B'))' = (OUT')' = OUT
 *
 * Bear in mind that (CLMUL(A<<1,B))' != (CLMUL(A,B)<<1)'
 *
 * Steps are:
 *
 * 1. Compute CLMUL(A, B) << 1
 * 2. Reduce with Q' (treated as 64bit) following the Doc. B approach
 * @param a __m512i internally is [A3:A2:A1:A0]
 * @param b __m512i internally is [B3:B2:B1:B0]
 * @return
 */
__m512i gfmul_reflected_avx512_parallel(
    const __m512i a, const __m512i b)
{
    /* algorithm by pH-Valiu based on Intel Doc. B Table 9 */
    // The goal is to perform gfmul_reflected on all 4 quad words (__m128i) at the same time
    // Therefore, each step must be seen as happening independept on each quad word



    // Step 1.0: Multiplication
    __m512i a0b0, a1b0, a0b1, a1b1;
    a0b0 = _mm512_clmulepi64_epi128(a, b, 0x00);
    a1b0 = _mm512_clmulepi64_epi128(a, b, 0x01);
    a0b1 = _mm512_clmulepi64_epi128(a, b, 0x10);
    a1b1 = _mm512_clmulepi64_epi128(a, b, 0x11);


    __m512i lo, hi, mid;
    // compute mid
    mid = _mm512_xor_si512(a1b0, a0b1);

    // compute high and low (for each quad word respectively)
    lo = _mm512_xor_si512(a0b0, _mm512_bslli_epi128(mid, 8));
    hi = _mm512_xor_si512(a1b1, _mm512_bsrli_epi128(mid, 8));


    // Step 1.1: Bitshift << 1
    __m512i lo_32bitMSB, hi_32bitMSB, lo_128bitMSB;
    lo_32bitMSB = _mm512_srli_epi32(lo, 31);
    hi_32bitMSB = _mm512_srli_epi32(hi, 31);
    lo = _mm512_slli_epi32(lo, 1);
    hi = _mm512_slli_epi32(hi, 1);
    lo_128bitMSB = _mm512_bsrli_epi128(lo_32bitMSB, 12);
    lo_32bitMSB = _mm512_bslli_epi128(lo_32bitMSB, 4);
    hi_32bitMSB = _mm512_bslli_epi128(hi_32bitMSB, 4);
    lo = _mm512_or_si512(lo, lo_32bitMSB);
    hi = _mm512_or_si512(hi, hi_32bitMSB);
    hi = _mm512_or_si512(hi, lo_128bitMSB);


    // Step 2.1: Reduce
    __m512i x = _mm512_clmulepi64_epi128(lo, Q_r_512, 0x00);
    hi = _mm512_xor_si512(hi, _mm512_bsrli_epi128(x, 8));
    hi = _mm512_xor_si512(hi, _mm512_unpacklo_epi64(lo, ZERO_512));

    lo = _mm512_xor_si512(lo, _mm512_bslli_epi128(x, 8));

    // Step 2.2: Reduce
    x = _mm512_clmulepi64_epi128(lo, Q_r_512, 0x01);
    hi = _mm512_xor_si512(hi, x);
    hi = _mm512_xor_si512(hi, _mm512_unpackhi_epi64(ZERO_512, lo));


    return hi;
}

/**
 * @brief print512_hex_lanes
 * This method prints them out in a way such that when you copy past it, to reuse the value, you may only reuse those values in a _mm512_set_epi32 instantiation or comparable.
 * Whereas if you want to load the bytes into a QByteArray, you have to apply BSWAP_MASK first
 * e.g.
 * print512_hex_lanes(var) -> "857119A1 F93A08AF D25E753D F3061F4B ... "
 *
 * Then to make a QByteArray out of it, you need to state: QByteArray::fromHex("... 4b1f06f33d755ed2af083af9a1197185");
 * You see how every byte is switched (LE<>BE)
 *
 * @param var
 * @return
 */
QString print512_hex_lanes(const __m512i& var)
{
    uint32_t val[16];
    memcpy(val, &var, sizeof(val));

    // print high lane first
    QString s = QString::asprintf(
        "%08X %08X %08X %08X %08X %08X %08X %08X %08X %08X %08X %08X %08X %08X %08X %08X",
        val[15], val[14], val[13], val[12], val[11], val[10], val[9], val[8], val[7], val[6], val[5], val[4], val[3], val[2], val[1], val[0]);
    return s;
}

bool isAligned64(const void* p) {
    return (reinterpret_cast<uintptr_t>(p) & 63) == 0;
}

void gfmul_reflected_avx512_parallel_test(){
    QByteArray a3_b = QByteArray::fromHex("9ca1a8db76a56facaea9510ce2e91845");
    QByteArray a2_b = QByteArray::fromHex("1d98752cf482a50136339f5e32f58d11");
    QByteArray a1_b = QByteArray::fromHex("9171b78cafdaed1399e4570c16b0e3a2");
    QByteArray a0_b = QByteArray::fromHex("50e27ae11262ca8c81f992482a947930");
    QByteArray b3_b = QByteArray::fromHex("2fd845fa7a215674f6de500b99dce012");
    QByteArray b2_b = QByteArray::fromHex("1a64fa8afa41d604a273670c413a7a39");
    QByteArray b1_b = QByteArray::fromHex("fc2962221dc1f7f0f5227b54c850ef96");
    QByteArray b0_b = QByteArray::fromHex("8c8a3eb3d5c7400abfbcee28678335c0");
    __m128i A3 = _mm_loadu_si128((__m128i*) a3_b.constData());
    __m128i A2 = _mm_loadu_si128((__m128i*) a2_b.constData());
    __m128i A1 = _mm_loadu_si128((__m128i*) a1_b.constData());
    __m128i A0 = _mm_loadu_si128((__m128i*) a0_b.constData());
    __m128i B3 = _mm_loadu_si128((__m128i*) b3_b.constData());
    __m128i B2 = _mm_loadu_si128((__m128i*) b2_b.constData());
    __m128i B1 = _mm_loadu_si128((__m128i*) b1_b.constData());
    __m128i B0 = _mm_loadu_si128((__m128i*) b0_b.constData());

    QByteArray a_b = QByteArray::fromHex("50e27ae11262ca8c81f992482a9479309171b78cafdaed1399e4570c16b0e3a21d98752cf482a50136339f5e32f58d119ca1a8db76a56facaea9510ce2e91845");
    QByteArray b_b = QByteArray::fromHex("8c8a3eb3d5c7400abfbcee28678335c0fc2962221dc1f7f0f5227b54c850ef961a64fa8afa41d604a273670c413a7a392fd845fa7a215674f6de500b99dce012");
    alignas(64) unsigned char ABytes[64];
    alignas(64) unsigned char BBytes[64];
    memcpy(ABytes, a_b.constData(), 64);
    memcpy(BBytes, b_b.constData(), 64);
    alignas(64) __m512i A = _mm512_load_si512((__m512i*) ABytes);
    alignas(64) __m512i B = _mm512_load_si512((__m512i*) BBytes);

    __m128i A3B3 = gfmul_reflected(A3, B3);
    __m128i A2B2 = gfmul_reflected(A2, B2);
    __m128i A1B1 = gfmul_reflected(A1, B1);
    __m128i A0B0 = gfmul_reflected(A0, B0);

    alignas(64) __m512i AB;
    gfmul_reflected_avx512(&A, &B, &AB);
    alignas(64) char t[64];
    _mm512_storeu_si512((__m512i*)t, AB);
    QByteArray ab_byteArray(t, 64);

    _mm_storeu_si128(&((__m128i*)t)[0], A0B0);
    _mm_storeu_si128(&((__m128i*)t)[1], A1B1);
    _mm_storeu_si128(&((__m128i*)t)[2], A2B2);
    _mm_storeu_si128(&((__m128i*)t)[3], A3B3);
    QByteArray ab_constructed_byteArray(t, 64);


    if(ab_byteArray == ab_constructed_byteArray){
        qInfo() << "[TEST - GFMUL AVX512] Assertion:"<< "OK";
    } else{
        qInfo() << "[TEST - GFMUL AVX512] Assertion: "<< "WRONG";
        qInfo() << "[TEST - GFMUL AVX512] Actual (normal): "<<ab_constructed_byteArray.toHex();
        qInfo() << "[TEST - GFMUL AVX512] Actual (512): "<<ab_byteArray.toHex();
    }

}
#endif


/**
 * @brief print128_hex_lanes
 * This method prints them out in a way such that when you copy past it, to reuse the value, you may only reuse those values in a _mm_set_epi32 instantiation or comparable.
 * Whereas if you want to load the bytes into a QByteArray, you have to apply BSWAP_MASK first
 * e.g.
 * print128_hex_lanes(var) -> "857119A1 F93A08AF D25E753D F3061F4B"
 *
 * Then to make a QByteArray out of it, you need to state: QByteArray::fromHex("4b1f06f33d755ed2af083af9a1197185");
 * You see how every byte is switched (LE<>BE)
 *
 * @param var
 * @return
 */
QString print128_hex_lanes(__m128i var)
{
    uint32_t val[4];
    memcpy(val, &var, sizeof(val));

    // print high lane first
    QString s = QString::asprintf(
        "%08X %08X %08X %08X",
        val[3], val[2], val[1], val[0]);
    return s;
}

#include <smmintrin.h>

#include <emmintrin.h>
#include <smmintrin.h>
#include <stdint.h>
#include <QDebug>


void gfmul_times_four_test(){
    // We are testing gfmul_times_four_reflected to 4 times gfmul_reflected
    __m128i H = _mm_loadu_si128((__m128i*) QByteArray::fromHex("b83b533708bf535d0aa6e52980d53b78").constData());
    __m128i H2 = gfmul_reflected(H, H);
    __m128i H3 = gfmul_reflected(H2, H);
    __m128i H4 = gfmul_reflected(H3, H);
    __m128i X = ZERO;
    __m128i X1 = _mm_loadu_si128((__m128i*) QByteArray::fromHex("2263f048f9ee49f51b22ff863a9e3b1f").constData());
    __m128i X2 = _mm_loadu_si128((__m128i*) QByteArray::fromHex("81a41759f340d74d82d270f084b8f522").constData());
    __m128i X3 = _mm_loadu_si128((__m128i*) QByteArray::fromHex("8967e442a24719c0a336ad05a025de40").constData());
    __m128i X4 = _mm_loadu_si128((__m128i*) QByteArray::fromHex("8998b973d729daf0beffbf01a88e5d5e").constData());


    // normal approach
    __m128i X_final = _mm_xor_si128(X, X1);
    X_final = gfmul_reflected(X_final, H);
    X_final = _mm_xor_si128(X_final, X2);
    X_final = gfmul_reflected(X_final, H);
    X_final = _mm_xor_si128(X_final, X3);
    X_final = gfmul_reflected(X_final, H);
    X_final = _mm_xor_si128(X_final, X4);
    X_final = gfmul_reflected(X_final, H);
    char t[16];
    _mm_storeu_si128((__m128i*)t, X_final);
    QByteArray X_final_normal_bytes(t, 16);

    // fast approach
    __m128i X_final_fast = gfmul_times_four_reflected(X4, X3, X2, X1, H, H2, H3, H4);
    _mm_storeu_si128((__m128i*)t, X_final_fast);
    QByteArray X_final_fast_bytes(t, 16);

    QByteArray assert_bytes = QByteArray::fromHex("4b1f06f33d755ed2af083af9a1197185");
    if(X_final_normal_bytes == assert_bytes && X_final_fast_bytes == assert_bytes){
        qInfo() << "[TEST - GFMUL TIMES FOUR] Assertion:"<< "OK";
    } else{
        qInfo() << "[TEST - GFMUL TIMES FOUR] Assertion: "<< "WRONG";
        qInfo() << "[TEST - GFMUL TIMES FOUR] Expected: "<<assert_bytes.toHex();
        qInfo() << "[TEST - GFMUL TIMES FOUR] Actual (normal): "<<X_final_normal_bytes.toHex();
        qInfo() << "[TEST - GFMUL TIMES FOUR] Actual (fast): "<<X_final_fast_bytes.toHex();
    }
}

void gfmul_test(){
    /* Test vektoren (TEST 1 - nach Intel Doc 2014 - Page 78 in Doc A.: "Intel Carry-Less Multiplication Instruction and its Usage for Computing in GCM mode"):
    / a: 7b5b5465 73745665 63746f72 5d53475d
    / b: 48692853 68617929 5b477565 726f6e5d
    / q: 00000000 00000000 00000000 00000087
    / assert: GFMUL(a,b,q) = c = 040229a0 9a5ed12e 7e4e10da 323506d2
    */
    __m128i a = _mm_set_epi32(0x7b5b5465, 0x73745665, 0x63746f72, 0x5d53475d);
    __m128i b = _mm_set_epi32(0x48692853, 0x68617929, 0x5b477565, 0x726f6e5d);
    __m128i res_assert = _mm_set_epi32(0x040229a0, 0x9a5ed12e, 0x7e4e10da, 0x323506d2);
    qInfo() << "a: "<<print128_hex_lanes(a)<<", b: "<<print128_hex_lanes(b);
    __m128i res = gfmul_k_optimized(a,b);
    qInfo() << "res: (a, b, q):\n|>"<<print128_hex_lanes(res);
    if(_mm_test_all_zeros(_mm_set_epi32(0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff), _mm_xor_si128(res_assert, res))) {
        qInfo() << "Assertion (res == res_assert): holds true!";
    } else {
        qWarning() << "Assertion (res == res_assert): is false!";
    }
    qInfo("\n");


    /* Test vektoren (TEST 2 - nach Intel Doc 2014 - Page 78 "...")
    / a: 952b2a56 a5604ac0 b32b6656 a05b40b6
    / b: dfa6bf4d ed81db03 ffcaff95 f830f061
    / q: 00000000 00000000 00000000 00000087
    / a_refl: 6d02da056a66d4cd035206a56a54d4a9
    / b_refl: 860f0c1fa9ff53ffc0db81b7b2fd65fb
    / assert: GFMUL(a_refl,b_refl,q) = c = 065B7FC3 340123F2 6DDAA34B 50D7CA5B      (actually, following our notation this would be c_refl)
    / assert: c_refl = da53eb0ad2c55bb64fc4802cc3feda60                             (and this would be c. But we followed the Intel naming scheme)
    / => In order for this to work, q must be q and must not be q_refl (this could also be due to us using a multiply&reduce routine different from Doc A.)
    */
    a = _mm_set_epi32(0x952b2a56, 0xa5604ac0, 0xb32b6656, 0xa05b40b6);
    b = _mm_set_epi32(0xdfa6bf4d, 0xed81db03, 0xffcaff95, 0xf830f061);
    __m128i a_refl = reflect_xmm(a);
    __m128i b_refl = reflect_xmm(b);
    __m128i c_assert = _mm_set_epi32(0x065B7FC3, 0x340123F2, 0x6DDAA34B, 0x50D7CA5B);
    __m128i c_refl_assert = _mm_set_epi32(0xda53eb0a, 0xd2c55bb6, 0x4fc4802c, 0xc3feda60);
    qInfo() << "a: "<<print128_hex_lanes(a)<<", b: "<<print128_hex_lanes(b);
    qInfo() << "a_refl: "<<print128_hex_lanes(a_refl)<<", b_refl: "<<print128_hex_lanes(b_refl);
    //qInfo() << "q_refl: "<<print128_hex_lanes(q_refl);
    __m128i c = gfmul_k_optimized(a_refl, b_refl);       // <<< switch here between gfmul and gfmul_k_optimized
    __m128i c_refl = reflect_xmm(c);

    qInfo() << "c: (a_refl, b_refl, q): \n|>"<<print128_hex_lanes(c);
    qInfo() << "c_refl: \n|>"<<print128_hex_lanes(c_refl);
    if(_mm_test_all_zeros(_mm_set_epi32(0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff), _mm_xor_si128(c_assert, c))) {
        qInfo() << "Assertion (c == c_assert): holds true!";
    } else {
        qWarning() << "Assertion (c == c_assert): is false!";
    }
    if(_mm_test_all_zeros(_mm_set_epi32(0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff), _mm_xor_si128(c_refl_assert, c_refl))) {
        qInfo() << "Assertion (c_refl == c_refl_assert): holds true!";
    } else {
        qWarning() << "Assertion (c_refl == c_refl_assert): is false!";
    }


    /*
     * Test of gfmul_reversed
     *
     * We try to get around using reflect_xmm, thats why OUT should already be the final value
     *
     * Assert against res: da53eb0a d2c55bb6 4fc4802c c3feda60
     * res_refl is: 065B7FC3 340123F2 6DDAA34B 50D7CA5B
     *
     */
    qInfo() <<"\nNow testing gfmul_reversed_bl_opt";
    a = _mm_set_epi32(0x952b2a56, 0xa5604ac0, 0xb32b6656, 0xa05b40b6);
    b = _mm_set_epi32(0xdfa6bf4d, 0xed81db03, 0xffcaff95, 0xf830f061);
    res_assert = _mm_set_epi32(0xda53eb0a, 0xd2c55bb6, 0x4fc4802c, 0xc3feda60);
    __m128i res_assert_refl = _mm_set_epi32(0x065B7FC3, 0x340123F2, 0x6DDAA34B, 0x50D7CA5B);
    qInfo() << "a: "<<print128_hex_lanes(a)<<", b: "<<print128_hex_lanes(b);
    gfmul_reflected_avx128(&a, &b, &res);
    //res = gfmul_reflected(a,b);
    __m128i res_refl = reflect_xmm(res);
    qInfo() << "res: (a, b, q):\n|>"<<print128_hex_lanes(res);
    //qInfo() << "res_refl: \n|>"<<print128_hex_lanes(res_refl);
    qInfo() << "res_assert:\n|>"<<print128_hex_lanes(res_assert);
    //qInfo() << "res_asser_refl:\n|>"<<print128_hex_lanes(res_assert_refl);
    if(_mm_test_all_zeros(_mm_set_epi32(0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff), _mm_xor_si128(res_assert, res))) {
        qInfo() << "Assertion (res == res_assert): holds true!";
    } else {
        qWarning() << "Assertion (res == res_assert): is false!";
    }
    qInfo("\n");


    uint8_t arr[16];
    uint8_t iv[12] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
    __m128i BSWAP_MASK = _mm_set_epi8(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);
    QByteArray d = QByteArray::fromRawData((char*) iv, 12);
    __m128i test2 = _mm_loadu_si128((__m128i*)d.constData());
    test2 = _mm_insert_epi32(test2, 0x01000000, 3);
    test2 = _mm_shuffle_epi8(test2, BSWAP_MASK);
    memcpy(arr, &test2, 16);
    qInfo() << "test: a15:|" <<arr[15]<<", "<<arr[14]<<", "<<arr[13]<<", "<<arr[12]<<", "<<arr[11]<<", "<<arr[10]<<", "<<arr[9]<<", "<<arr[8]<<", "<<arr[7]<<", "<<arr[6]<<", "<<arr[5]<<", "<<arr[4]<<", "<<arr[3]<<", "<<arr[2]<<", "<<arr[1]<<", "<<arr[0]<<"|:a0";

    __m128i test3 = _mm_set_epi32(0x01000000, *(const int*)(d.constData() + 8), *(const int*)(d.constData() + 4), *(const int*)(d.constData() + 0));
    memcpy(arr, &test3, 16);
    qInfo() << "test: a15:|" <<arr[15]<<", "<<arr[14]<<", "<<arr[13]<<", "<<arr[12]<<", "<<arr[11]<<", "<<arr[10]<<", "<<arr[9]<<", "<<arr[8]<<", "<<arr[7]<<", "<<arr[6]<<", "<<arr[5]<<", "<<arr[4]<<", "<<arr[3]<<", "<<arr[2]<<", "<<arr[1]<<", "<<arr[0]<<"|:a0";


    const __m128i x = _mm_set_epi32(0x952b2a56, 0xa5604ac0, 0xb32b6656, 0xa05b40b6);
    const __m128i y = _mm_set_epi32(0xdfa6bf4d, 0xed81db03, 0xffcaff95, 0xf830f061);
    //gfmul_reflected_avx512_parallel(x512, x512);
    gfmul_times_four_test();
    gfmul_reflected_avx512_parallel_test();

    __m512i x512 = _mm512_set_epi32(0x952b2a56, 0x952b2a56, 0x952b2a56, 0x952b2a56, 0x952b2a56, 0x952b2a56, 0x952b2a56, 0x952b2a56, 0x952b2a56, 0x952b2a56, 0x952b2a56, 0x952b2a56, 0x952b2a56, 0x952b2a56, 0x952b2a56, 0x952b2a56);
    __m512i result;
    gfmul_reflected_avx512(&x512, &x512, &result);
    BenchmarkUtil::run("gfmul_reflected", gfmul_reflected, x, y);
}
