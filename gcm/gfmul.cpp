#include <gcm/gfmul.h>
#include "gcm/benchmarkutil.h"
#include <QDebug>


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

__m128i gfmul(__m128i a, __m128i b){
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
__m128i gfmul_k_optimized(__m128i a, __m128i b){
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
__m128i gfmul_reflected(__m128i a, __m128i b){
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

int fun() {
    uint8_t iv[12] = { 1, 2, 3, 4,  5, 6, 7, 8,  9, 10, 11, 12 };

    // --- CASE A: treat each IV chunk as LITTLE-ENDIAN 32-bit integers ---
    __m128i a = _mm_set_epi32(
        0x01000000,
        *(const int*)(iv + 8),
        *(const int*)(iv + 4),
        *(const int*)(iv + 0)
        );

    uint8_t outA[16];
    memcpy(outA, &a, 16);

    qInfo() << "\n=== CASE A (little-endian int loads) ===";
    qInfo() << QByteArray((char*)outA,16).toHex(' ');


    // --- CASE B: treat IV as BIG-ENDIAN byte string (GCM correct) ---
    // Load bytes directly in the correct order
    uint8_t block[16] = {
        iv[0],iv[1],iv[2],iv[3],
        iv[4],iv[5],iv[6],iv[7],
        iv[8],iv[9],iv[10],iv[11],
        0x00,0x00,0x00,0x01  // 32-bit counter in big-endian
    };

    __m128i b = _mm_loadu_si128((__m128i*)block);

    uint8_t outB[16];
    memcpy(outB, &b, 16);

    qInfo() << "\n=== CASE B (byte-wise big-endian load, GCM correct) ===";
    qInfo() << QByteArray((char*)outB,16).toHex(' ');

    return 0;
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
    res = gfmul_reflected(a,b);
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

    fun();
}
