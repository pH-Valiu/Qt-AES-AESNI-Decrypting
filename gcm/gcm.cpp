#include <gcm/gcm.h>
#include <QDebug>
#include <QString>

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

__m128i gfmul(__m128i a, __m128i b, __m128i q){
    // Step 1: Multiply
    __m128i a0b0 = _mm_clmulepi64_si128(a, b, 0x00);
    __m128i a0b1 = _mm_clmulepi64_si128(a, b, 0x10);
    __m128i a1b0 = _mm_clmulepi64_si128(a, b, 0x01);
    __m128i a1b1 = _mm_clmulepi64_si128(a, b, 0x11);

    __m128i mid = _mm_xor_si128(a0b1, a1b0); // computes mid = A0B1 + A1B0

    __m128i c01 = _mm_xor_si128(a0b0, _mm_slli_si128(mid, 8)); // computes C[1:0] = A0B0 + (mid << x^64)
    __m128i c23 = _mm_xor_si128(a1b1, _mm_srli_si128(mid, 8)); // computes C[3:2] = A1B1 + (mid >> x^64)

    // Step 2.1: Reduce
    __m128i x = _mm_clmulepi64_si128(a1b1, q, 0x01);  // computes C[3] * Q = upper(A1B1) * lower(Q)
    c01 = _mm_xor_si128(c01, _mm_slli_si128(x, 8)); // add lower half of x (X[0]) to upper part of C[1:0]
    c23 = _mm_xor_si128(c23, _mm_srli_si128(x, 8)); // add higher half of x (X[1]) to lower part of C[3:2] (higher part is just dangling around)

    // Step 2.2: Reduce
    x = _mm_clmulepi64_si128(c23, q, 0x00);         // works because higher part is not used in the calculation
    c01 = _mm_xor_si128(c01, x);

    return c01;
}

__m128i gfmul(__m128i a, __m128i b){
    __m128i q = _mm_set_epi32(0, 0, 0, 0x00000087); //0x87 is 10000111 being x^7 + x^2 + x + 1

    return gfmul(a,b, q);
}

__m128i bitshift_left(__m128i a, unsigned char count){
    __m128i carry = _mm_slli_si128(a, 8);   // old compilers only have the confusingly named _mm_slli_si128 synonym
    if (count >= 64)
        return _mm_slli_epi64(carry, count-64);  // the non-carry part is all zero, so return early
    // else
    carry = _mm_srli_epi64(carry, 64-count);  // After bslli shifted left by 64b

    a = _mm_slli_epi64(a, count);
    return _mm_or_si128(a, carry);
}

#include <cpuid.h>
#include <iostream>

bool hasSSSE3() {
    unsigned int eax, ebx, ecx, edx;
    if (__get_cpuid(1, &eax, &ebx, &ecx, &edx)) {
        return (ecx & bit_SSSE3);
    }
    return false;
}

QString print128_hex_lanes(__m128i var)
{
    uint32_t val[4];
    memcpy(val, &var, sizeof(val));

    // print high lane first
    QString s = QString::asprintf(
        "%08X %08X %08X %08X\n",
        val[3], val[2], val[1], val[0]);
    return s;
}

void gfmul_test(){
    qInfo() << "hasSSSE3?: " << hasSSSE3();
    /* Test vektoren (TEST 1 - nach Intel Doc 2014 - Page 78 "Intel Carry-Less Multiplication Instruction and its Usage for Computing in GCM mode"):
    / a: 7b5b5465 73745665 63746f72 5d53475d
    / b: 48692853 68617929 5b477565 726f6e5d
    / q: 00000000 00000000 00000000 00000087
    / assert: GFMUL(a,b,q) = 040229a0 9a5ed12e 7e4e10da 323506d2
    */
    __m128i a = _mm_set_epi32(0x7b5b5465, 0x73745665, 0x63746f72, 0x5d53475d);
    __m128i b = _mm_set_epi32(0x48692853, 0x68617929, 0x5b477565, 0x726f6e5d);
    __m128i q = _mm_set_epi32(0x00000000, 0x00000000, 0x00000000, 0x00000087);
    qInfo() << "a: "<<print128_hex_lanes(a)<<", b: "<<print128_hex_lanes(b);
    __m128i res = gfmul(a,b,q);
    qInfo() << "res: (a, b, q):\n|>"<<print128_hex_lanes(res)<<"\n";


    /* Test vektoren (TEST 2 - nach Intel Doc 2014 - Page 78 "...")
    / a: 952b2a56 a5604ac0 b32b6656 a05b40b6
    / b: dfa6bf4d ed81db03 ffcaff95 f830f061
    / q: 00000000 00000000 00000000 00000087
    / a_refl: 6d02da056a66d4cd035206a56a54d4a9
    / b_refl: 860f0c1fa9ff53ffc0db81b7b2fd65fb
    / assert: GFMUL(a_refl,b_refl,q) = c = 65B7FC3340123F26DDAA34B50D7CA5B
    / assert: c_refl = da53eb0ad2c55bb64fc4802cc3feda60
    / => In order for this to work, q must be q and must not be q_refl
    */
    a = _mm_set_epi32(0x952b2a56, 0xa5604ac0, 0xb32b6656, 0xa05b40b6);
    b = _mm_set_epi32(0xdfa6bf4d, 0xed81db03, 0xffcaff95, 0xf830f061);
    q = _mm_set_epi32(0x00000000, 0x00000000, 0x00000000, 0x00000087);
    __m128i a_refl = reflect_xmm(a);
    __m128i b_refl = reflect_xmm(b);
    qInfo() << "a: "<<print128_hex_lanes(a)<<", b: "<<print128_hex_lanes(b);
    qInfo() << "a_refl: "<<print128_hex_lanes(a_refl)<<", b_refl: "<<print128_hex_lanes(b_refl);
    //qInfo() << "q_refl: "<<print128_hex_lanes(q_refl);
    __m128i c = gfmul(a_refl, b_refl, q);
    __m128i c_refl = reflect_xmm(c);

    qInfo() << "c: (a_refl, b_refl, q): \n|>"<<print128_hex_lanes(c);
    qInfo() << "c_refl: \n|>"<<print128_hex_lanes(c_refl);

}



void encrypt(){
    /*
     * 1. take params:
     *  - message, key, iv, (aad)
     *  ---> output: ciphertext, tag
     *
     * 2. key expansions for all 10/11/12 AES key rounds
     *
     * 4. compute Y0 consisting of nonce (iv with 96bit or ghash(iv)) and 0x00000001),
     *  - 4.1 if IV not 96 bits, let s be the number of bits missing to append to IV to make it to full 128bits blocks. s =128 * upper(len(IV)/128) - len(IV),
     *  -     then Y0 = GHASH_H(IV||0^{s+64}||[len(IV)]_64) (actually: multiples of the block size)
     *  - 4.2 IV construction is handled by two methods: (see NIST 800-38d) is supported IV lengths < 96, then deterministic approach with a device distinct fixed field (MSB-32bits), and an invocation field(LSB-64bits), being an integer counter of linear feedback shift register with adequate cycling.
     *  -     for supported IV lengths >= 96, RBG based construction shall be used with random field (>= 96bits) and a free free field
     *  - 4.3 IV construction methods must be supplied
     *  - 4.4 side note: no more than 2^32 with the same key (given by IV would be reused due to 32bit 0x00000001 field)
     *  - 4.5 side note: a loss of power shall not cause the repetition of the IVs. (could be realised by storing the IV invocation field in non-volatile memory already a few states ahead)
     *  - 4.6 side note: in NIST 800-38d, Appendix B there are warning issued for using AES-GCM with long message blocks as each additional block increases to probability for a targeted ciphertext forgery to succeed. --> adhere to IV, message, key size limits!!!
     * 5. compute H being 0x0000... encrypted with AES and key
     * 6. compute keystream using Y0+i for as many bits required as message is long (AES-CTR mode needs no padding). The first Y0+i used, is Y1, corresponding to an Y value with a 2 at the end.
     * 7. create ciphertext using XOR
     * 8. compute tag Xi using X0=0x0 and Xi=(Xi-1 XOR Ci) * H  (later we will refine it by doing blockwise x8 or x4 approach. requires precomputation of H^2, H^3,...)
     *  - 8.1, the full input to the tag generation is (A||0^v||C||0^u||[len(A)]_64||[len(C)]_64) with v and u being the respective number of bits missing A and C to fill complete 128bits blocks. (actually: multiples of the block size)
     *  - 8.2 the final Tag is then T = MSB_t(GCTR_K(Y0, S)), where t is tag size, and S is the final value of the GHASH computation scheme.
     *
     * side note: GHASH specifies all operandi to be bit-reflected
     *
     */
}
