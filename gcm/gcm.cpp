#include <gcm/gcm.h>
#include <QDebug>
#include <QString>

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

struct int256 reflect_ymm(struct int256 X){
    struct int256 R;

    R.t10 = reflect_xmm(X.t32);
    R.t32 = reflect_xmm(X.t10);

    return R;
}

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
__m128i gfmul_k_optimized(__m128i a, __m128i b, __m128i q){
    // Step 0: Pre-requesites
    __m128i k = _mm_xor_si128(gfmul(_mm_srli_si128(b, 8), q, q), _mm_slli_si128(b, 8));     // K = GFMUL(B[1], Q) + B[0]*x^64

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
    __m128i x = _mm_clmulepi64_si128(c23, q, 0x00);
    c01 = _mm_xor_si128(x, c01);

    return c01;
}

__m128i gfmul_k_optimized(__m128i a, __m128i b){
    __m128i q = _mm_set_epi32(0, 0, 0, 0x00000087); //0x87 is 10000111 being x^7 + x^2 + x + 1
    return gfmul_k_optimized(a, b, q);
}

__m128i gfmul_k_optimized_reversed(__m128i a, __m128i b, __m128i q){
    __m128i a_refl = reflect_xmm(a);
    __m128i b_refl = reflect_xmm(b);
    //__m128i q_refl = reflect_xmm(bitshift_right(q, 1));
    __m128i q_refl = _mm_set_epi32(0, 0, 0xc2000000, 0);        // this is reflect_64(Q >> 1)
    __m128i k = _mm_xor_si128(gfmul(_mm_srli_si128(b, 8), q, q), _mm_slli_si128(b, 8));     // K = GFMUL(B[1], Q) + B[0]*x^64
    __m128i k_refl = reflect_xmm(k);

    // Step 1: Multiply
    __m128i a0k0 = _mm_clmulepi64_si128(a_refl, k_refl, 0x00);
    __m128i a0k1 = _mm_clmulepi64_si128(a_refl, k_refl, 0x10);
    __m128i a1b0 = _mm_clmulepi64_si128(a_refl, b_refl, 0x01);
    __m128i a1b1 = _mm_clmulepi64_si128(a_refl, b_refl, 0x11);

    __m128i lower  = _mm_xor_si128(a1b0, a0k0);     // computes lower  = a1b0 + a0k0
    __m128i higher = _mm_xor_si128(a1b1, a0k1);     // computes higher = a1b1 + a0k1

    __m128i c01 = _mm_xor_si128(lower, _mm_slli_si128(higher, 8));  // c01 = lower + (higher << 64)
    __m128i c23 = _mm_srli_si128(higher, 8);                        // c02 = higher >> 64           (highest 64 bits are 0)

    // Step 2: Reduce
    __m128i x = _mm_clmulepi64_si128(c01, q_refl, 0x00);        // this is x = CLMUL(c0, Q')

    __m128i f = _mm_srli_si128(c01, 8);             // f = [0, c1]
    f = _mm_xor_si128(f, x);                        // f = [x1, c1+x0]
    f = _mm_xor_si128(f, _mm_slli_si128(c23, 8));   // f = [x1+c2, c1+x0]
    f = _mm_xor_si128(f, _mm_slli_si128(c01, 8));    // f = [x1+c2+c0, c1+x0]

    return f;
}

__m128i gfmul_k_optimized_reversed(__m128i a, __m128i b){
    __m128i q = _mm_set_epi32(0, 0, 0, 0x00000087); //0x87 is 10000111 being x^7 + x^2 + x + 1
    return gfmul_k_optimized_reversed(a, b, q);
}

struct int256 clmul(__m128i a, __m128i b){
    // Step 1: Multiply
    __m128i a0b0 = _mm_clmulepi64_si128(a, b, 0x00);
    __m128i a0b1 = _mm_clmulepi64_si128(a, b, 0x10);
    __m128i a1b0 = _mm_clmulepi64_si128(a, b, 0x01);
    __m128i a1b1 = _mm_clmulepi64_si128(a, b, 0x11);

    __m128i mid = _mm_xor_si128(a0b1, a1b0); // computes mid = A0B1 + A1B0

    __m128i c01 = _mm_xor_si128(a0b0, _mm_slli_si128(mid, 8)); // computes C[1:0] = A0B0 + (mid << x^64)
    __m128i c23 = _mm_xor_si128(a1b1, _mm_srli_si128(mid, 8)); // computes C[3:2] = A1B1 + (mid >> x^64)

    return {c01, c23};
}


/**
 * This approach mostly strictly follows the "Simple Bit-Reflected GHASH Algorithm" of Table 9 in Doc B.
 * Important note: while the table says that inputs a and b are reflected, they are indeed not.
 * Additionally, after the CLMUL(A,B), the result must be shifted by one to the left (CLMUL(A,B) << 1) to accomadate for the fact that we are intrinsically working on reflected values.
 * Meanwhile, the reduction polynom is correctly Q' = reflect_64(Q>>1) = epi32(0...0, 0...0, 0xc2000000, 0...0)
 * The returnted value must not be reflected but instead immediately represents the correct output as if we had done:
 * OUT = (CLMUL(A', B'))' = (OUT')' = OUT
 *
 * Bear in mind that (CLMUL(A<<1,B))' != (CLMUL(A,B)<<1)'
 *
 * Steps are:
 *
 * 1. Compute CLMUL(A, B) << 1
 * 2. Reduce with Q' (treated as 64bit) following the Doc. B approach
 *
 * @brief gfmul_idea
 * @param a
 * @param b
 * @return
 */
__m128i gfmul_reversed(__m128i a, __m128i b){
    __m128i q = _mm_set_epi32(0, 0, 0xc2000000, 0);

    // Step 1: Multiply
    __m128i a0b0 = _mm_clmulepi64_si128(a, b, 0x00);
    __m128i a0b1 = _mm_clmulepi64_si128(a, b, 0x10);
    __m128i a1b0 = _mm_clmulepi64_si128(a, b, 0x01);
    __m128i a1b1 = _mm_clmulepi64_si128(a, b, 0x11);

    __m128i mid = _mm_xor_si128(a0b1, a1b0);      // computes mid = A0B1 + A1B0

    __m128i c01 = _mm_xor_si128(a0b0, _mm_slli_si128(mid, 8));    // computes C[1:0] = A0B0 + (mid << x^64)
    __m128i c23 = _mm_xor_si128(a1b1, _mm_srli_si128(mid, 8));    // computes C[3:2] = A1B1 + (mid >> x^64)

    // Step 1.1: Bitshift << 1
    struct int256 bsl = bitshift_left256(c23, c01, 1);
    c01 = bsl.t10;
    c23 = bsl.t32;

    // Step 2.1: Reduce
    __m128i x = _mm_clmulepi64_si128(c01, q, 0x00);
    c23 = _mm_xor_si128(c23, _mm_srli_si128(x, 8));       // add higher half of x (X[1]) to lower part of C[3:2]
    c23 = _mm_xor_si128(c23, _mm_and_si128(c01, _mm_set_epi64x(0, -1))); // add only C[0] to lower part of C[3:2]  (zeroing out C[1])

    c01 = _mm_xor_si128(c01, _mm_slli_si128(x, 8));       // add lower half of x (X[0]) to upper part of C[1:0]

    // Step 2.2: Reduce
    x = _mm_clmulepi64_si128(c01, q, 0x01);               // computes C[1] * Q = higher(C[1:0]) * lower(Q)
    c23 = _mm_xor_si128(c23, x);                          // add full X on C[3:2]
    c23 = _mm_xor_si128(c23, _mm_and_si128(c01, _mm_set_epi64x(-1, 0)));  // add C[1] on higher C[3:2] (zeroing out C[0])

    return c23;
}

/**
 * @brief gfmul_reversed_bl_opt: This version is just the same as gfmul_reversed, just with the addition of the bitshift operation being optimized
 * instead of using our custom own bitshift_left and _right operations.
 * @param a
 * @param b
 * @return
 */
__m128i gfmul_reversed_bl_opt(__m128i a, __m128i b){
    //__m128i q = _mm_set_epi32(0, 0, 0xc2000000, 0);

    // Step 1: Multiply
    __m128i a0b0 = _mm_clmulepi64_si128(a, b, 0x00);
    __m128i a0b1 = _mm_clmulepi64_si128(a, b, 0x10);
    __m128i a1b0 = _mm_clmulepi64_si128(a, b, 0x01);
    __m128i a1b1 = _mm_clmulepi64_si128(a, b, 0x11);

    __m128i mid = _mm_xor_si128(a0b1, a1b0);      // computes mid = A0B1 + A1B0

    __m128i c01 = _mm_xor_si128(a0b0, _mm_slli_si128(mid, 8));    // computes C[1:0] = A0B0 + (mid << x^64)
    __m128i c23 = _mm_xor_si128(a1b1, _mm_srli_si128(mid, 8));    // computes C[3:2] = A1B1 + (mid >> x^64)

    // Step 1.1: Bitshift << 1
    __m128i tmp7,tmp8,tmp9;
    tmp7 = _mm_srli_epi32(c01, 31);
    tmp8 = _mm_srli_epi32(c23, 31);
    c01 = _mm_slli_epi32(c01, 1);
    c23 = _mm_slli_epi32(c23, 1);
    tmp9 = _mm_srli_si128(tmp7, 12);
    tmp8 = _mm_slli_si128(tmp8, 4);
    tmp7 = _mm_slli_si128(tmp7, 4);
    c01 = _mm_or_si128(c01, tmp7);
    c23 = _mm_or_si128(c23, tmp8);
    c23 = _mm_or_si128(c23, tmp9);


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

__m128i gfmul_reversed_k_optimized(__m128i a, __m128i b){
    /*
     * Scrabbling idea
     * I try to correct DocB gfmul k optimized.
     * The K definition is wrong I think. Because they simply use the reversed K' = reflect(K), but I dont think thats possible.
     * And I also dont want to reflect.
     * And I need to keep in mind this <<1 to the whole C[3:0].
     * But while in the end we do not have a C[3:0] but only a Y[2:0] left, we still have to think about how the K for our purposes is being created
     * K, in the normal world, is derived from seeing that the highest 64 bits = A[1] times B[1], when we want to reduce them, so:
     * GFMUL(A[1],B[1]) = GFMUL(GFMUL(A[1],B[1]), Q) mod P ~= A[1] * B[1] * Q mod P = A[1] * K = GFMUL(A[1], K)
     * here, we defined K = GFMUL(B[1], Q)
     * so: GFMUL(A[1], K) = CLMUL(A[1], K[0]) + CLMUL(A[1], K[1])*x^64
     *
     * mapping this into reversed world, we would have to look at GFMUL(A[0], B[0]) = GFMUL(GFMUL(A[0], B[0]), Q_r) mod P ?= GFMUL(A[0], K_r)
     * here, K_r = GFMUL(B[0], Q_r)
     * But (also in the top case), because B[0] (and B[1]), and Q_r (and Q) are all just 64bit big, Multiplying them together will never leave the finite field and thus require a reduction.
     * Hence, K_r = GFMUL(B[0], Q_r) == CLMUL(B[0], Q_r)
     * And because of the CL-identity, I think we would have to perform a shift: --> (CLMUL(B[0], Q_r) << 1). The carry out might have to be treated though.
     *
     * Additional food for though: The bitshift to B[0]. Must it be applied before CLMUL(B[0], Q_r) or after the CLMUL operation?
     * 1. We try after the CLMUL operation
     *
     *
     * Now, the second part will temporarily not be discussed, because one step at a time
     * Some idea though:
     * K_r = (CLMUL(B[0], Q_r) + B[1]*x^64) << 1 = (B[1]*x^64 + CLMUL(B[0], Q_r)) << 1
     * Still, there is the question of the carry-out
     * I would argue, that since we have remaining [A[1]*B[1]]*x^128 + [A[1]*B[0]]*x^64, we have to shift these M[2:0] << 1 and bring that carry out from before and set at LSB.
     *
     */
    __m128i q_k = Q_r;
    __m128i k_r = _mm_clmulepi64_si128(b, q_k, 0x00); // this is the same as GFMUL(B[0], Q_r) ?= CLMUL(B[0], Q_r)
    k_r = _mm_xor_si128(k_r, _mm_slli_si128(b, 8));     // this adds B[0]*x^64 onto K_r

    // we have to get the carry-out of k_r (when shifted by << 1)
    __m128i carry_out = _mm_srli_epi32(k_r, 31);
    carry_out = _mm_srli_si128(carry_out, 12);      // carry now contains only the most significant bit of k_r at index 0

    __m128i k_r_bl = bitshift_left(k_r, 1);     // now, we have shifted it once to the left (CLMUL(B[0], Q_r) << 1)

    // now, we would have to reduce k_r with q_k, but since 64bit times 64bit never leaves our F(2^128) field, we do not need to reduce



    // now that we have k_r, we compute the rest
    __m128i a0b1 = _mm_clmulepi64_si128(a, b, 0x10);
    __m128i a1b0 = _mm_clmulepi64_si128(a, b, 0x01);
    __m128i a1b1 = _mm_clmulepi64_si128(a, b, 0x11);
    __m128i a1k0 = _mm_clmulepi64_si128(a, k_r_bl, 0x01);       // A1K0 = CLMUL(A[1], K_r[0]);
    __m128i a1k1 = _mm_clmulepi64_si128(a, k_r_bl, 0x11);       // A1K1 = CLMUL(A[1], K_r[1]);

    __m128i mid = _mm_xor_si128(a0b1, a1b0);
    mid = _mm_xor_si128(mid, a1k0);

    __m128i c01 = _mm_slli_si128(mid, 8);
    __m128i c23 = _mm_xor_si128(a1b1, a1k1);
    c23 = _mm_xor_si128(c23, _mm_srli_si128(mid, 8));


    // bitshift << 1
    struct int256 t = bitshift_left256(c23, c01, 1);
    c01 = t.t10;
    c23 = t.t32;

    // Reduce with Q_r

    __m128i x = _mm_clmulepi64_si128(c01, q_k, 0x01);
    c23 = _mm_xor_si128(c23, x);
    c23 = _mm_xor_si128(c23, _mm_unpackhi_epi64(ZERO, c01));


    return c23;
}

__m128i gfmul_docA(__m128i a, __m128i b){

    // Step 1: Multiply
    __m128i a0b0 = _mm_clmulepi64_si128(a, b, 0x00);
    __m128i a0b1 = _mm_clmulepi64_si128(a, b, 0x10);
    __m128i a1b0 = _mm_clmulepi64_si128(a, b, 0x01);
    __m128i a1b1 = _mm_clmulepi64_si128(a, b, 0x11);

    __m128i mid = _mm_xor_si128(a0b1, a1b0);      // computes mid = A0B1 + A1B0

    __m128i c01 = _mm_xor_si128(a0b0, _mm_slli_si128(mid, 8));    // computes C[1:0] = A0B0 + (mid << x^64)
    __m128i c23 = _mm_xor_si128(a1b1, _mm_srli_si128(mid, 8));    // computes C[3:2] = A1B1 + (mid >> x^64)

    // Step 1.1: Shift << 1
    struct int256 bsl = bitshift_left256(c23, c01, 1);
    c01 = bsl.t10;
    c23 = bsl.t32;

    quint64 x0 = _mm_cvtsi128_si64(c01);
    quint64 x1 = _mm_cvtsi128_si64(_mm_srli_si128(c01, 8));
    quint64 x2 = _mm_cvtsi128_si64(c23);
    quint64 x3 = _mm_cvtsi128_si64(_mm_srli_si128(c23, 8));

    quint64 ai = x0 << 63;
    quint64 bi = x0 << 62;
    quint64 ci = x0 << 57;

    quint64 d = x1 ^ ai ^ bi ^ ci;

    __m128i dx0 = _mm_set_epi64x(d, x0);
    qInfo()<<"BREAK: [d = x1 ^ a ^ b ^ c, x0]\n>>>|"<<print128_hex_lanes(dx0);
    __m128i e10 = bitshift_right(dx0, 1);
    __m128i f10 = bitshift_right(dx0, 2);
    __m128i g10 = bitshift_right(dx0, 7);

    quint64 e0 = _mm_cvtsi128_si64(e10);
    quint64 e1 = _mm_cvtsi128_si64(_mm_srli_si128(e10, 8));
    quint64 f0 = _mm_cvtsi128_si64(f10);
    quint64 f1 = _mm_cvtsi128_si64(_mm_srli_si128(f10, 8));
    quint64 g0 = _mm_cvtsi128_si64(g10);
    quint64 g1 = _mm_cvtsi128_si64(_mm_srli_si128(g10, 8));

    c23 = _mm_xor_si128(c23, dx0);
    c23 = _mm_xor_si128(c23, e10);
    c23 = _mm_xor_si128(c23, f10);
    c23 = _mm_xor_si128(c23, g10);

    return c23;
}

__m128i gfmul_original_docA (__m128i a, __m128i b){
    __m128i tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9;

    // Test: instead of doing C[3:0] << 1, we could simply do A<<1 before CLMUL(A<<1,B)
    // --> Test result: NO! --> (CLMUL(A, B) << 1)' != (CLMUL(A<<1, B))'

    __m128i a0b0 = _mm_clmulepi64_si128(a, b, 0x00);    //3
    tmp3 = a0b0;
    __m128i a0b1 = _mm_clmulepi64_si128(a, b, 0x10);    //4
    tmp4 = a0b1;
    __m128i a1b0 = _mm_clmulepi64_si128(a, b, 0x01);    //5
    tmp5 = a1b0;
    __m128i a1b1 = _mm_clmulepi64_si128(a, b, 0x11);    //6
    tmp6 = a1b1;

    __m128i mid  = _mm_xor_si128(a0b1, a1b0);           //4
    tmp4 = mid;
    __m128i c01 = _mm_xor_si128(a0b0, _mm_slli_si128(mid, 8));    //3
    tmp3 = c01;
    __m128i c23 = _mm_xor_si128(a1b1, _mm_srli_si128(mid, 8));    //6
    tmp6 = c23;

    __m128i c01_shiftR31_epi32 = _mm_srli_epi32(tmp3, 31);      // carry-out bit of c01 of each 32bit lane
    tmp7 = c01_shiftR31_epi32;
    __m128i c23_shiftR31_epi32 = _mm_srli_epi32(tmp6, 31);      // carry-out bit of c23 of each 32bit lane
    tmp8 = c23_shiftR31_epi32;
    __m128i c01_shiftL1_epi32 = _mm_slli_epi32(tmp3, 1);        // shift <<1 of c01 in each 32 bit lane
    tmp3 = c01_shiftL1_epi32;
    __m128i c23_shiftL1_epi32 = _mm_slli_epi32(tmp6, 1);        // shift <<1 of c23 in each 32 bit lane
    tmp6 = c23_shiftL1_epi32;
    __m128i c01_shiftR31_epi32_shiftR96_si128 = _mm_srli_si128(tmp7, 12);   // contains just the carry of c01 (e.g. the highest bit, but now at position 0)
    tmp9 = c01_shiftR31_epi32_shiftR96_si128;
    __m128i c23_shiftR31_epi32_shiftL32_si128 = _mm_slli_si128(tmp8, 4);    // contains just the carries of c23 (32bit lane) but offset by one lane to the left: e.g. epi32(0...0carry_2, 0...0carry_1, 0...0carry_0, 0...0) of c23
    tmp8 = c23_shiftR31_epi32_shiftL32_si128;
    __m128i c01_shiftR31_epi32_shiftL32_si128 = _mm_slli_si128(tmp7, 4);    // contains just the carries of c01 (32bit lane) but offset by one lane to the left: e.g. epi32(0...0carry_2, 0...0carry_1, 0...0carry_0, 0...0) of c01
    tmp7 = c01_shiftR31_epi32_shiftL32_si128;

    __m128i c01_shiftL1_epi32_OR_c01_shiftR31_epi32_shiftL32_si128 = _mm_or_si128(tmp3, tmp7);      // comined all of c01
    tmp3 = c01_shiftL1_epi32_OR_c01_shiftR31_epi32_shiftL32_si128;
    __m128i c23_shiftL1_epi32_OR_c23_shiftR31_epi32_shiftL32_si128 = _mm_or_si128(tmp6, tmp8);
    tmp6 = c23_shiftL1_epi32_OR_c23_shiftR31_epi32_shiftL32_si128;
    __m128i c23_shiftL1_epi32_OR_c23_shiftR31_epi32_shiftL32_si128_OR_c01_shiftR31_epi32_shiftR96_si128 = _mm_or_si128(tmp6, tmp9);
    tmp6 = c23_shiftL1_epi32_OR_c23_shiftR31_epi32_shiftL32_si128_OR_c01_shiftR31_epi32_shiftR96_si128;

    // All the above does simply [tmp6:tmp3] << 1 = C[3:0] << 1
    // We have a command for that:
    struct int256 ret = bitshift_left256(c23, c01, 1);
    __m128i c23_sl1 = ret.t32;
    __m128i c01_sl1 = ret.t10;
    tmp6 = c23_sl1;
    tmp3 = c01_sl1;


    tmp7 = _mm_slli_epi32(tmp3, 31);
    tmp8 = _mm_slli_epi32(tmp3, 30);
    tmp9 = _mm_slli_epi32(tmp3, 25);
    tmp7 = _mm_xor_si128(tmp7, tmp8);
    tmp7 = _mm_xor_si128(tmp7, tmp9);
    tmp8 = _mm_srli_si128(tmp7, 4);
    tmp7 = _mm_slli_si128(tmp7, 12);
    tmp3 = _mm_xor_si128(tmp3, tmp7);
    tmp2 = _mm_srli_epi32(tmp3, 1);
    tmp4 = _mm_srli_epi32(tmp3, 2);
    tmp5 = _mm_srli_epi32(tmp3, 7);
    tmp2 = _mm_xor_si128(tmp2, tmp4);
    tmp2 = _mm_xor_si128(tmp2, tmp5);
    tmp2 = _mm_xor_si128(tmp2, tmp8);
    tmp3 = _mm_xor_si128(tmp3, tmp2);
    tmp6 = _mm_xor_si128(tmp6, tmp3);

    return tmp6;
}

__m128i gfmul_original_docA_mod(__m128i a, __m128i b){
    __m128i tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9;
    tmp3 = _mm_clmulepi64_si128(a, b, 0x00);
    tmp4 = _mm_clmulepi64_si128(a, b, 0x10);
    tmp5 = _mm_clmulepi64_si128(a, b, 0x01);
    tmp6 = _mm_clmulepi64_si128(a, b, 0x11);
    tmp4 = _mm_xor_si128(tmp4, tmp5);
    tmp5 = _mm_slli_si128(tmp4, 8);
    tmp4 = _mm_srli_si128(tmp4, 8);
    tmp3 = _mm_xor_si128(tmp3, tmp5);
    tmp6 = _mm_xor_si128(tmp6, tmp4);

    struct int256 bsl = bitshift_left256(tmp6, tmp3, 1);
    tmp6 = bsl.t32;
    tmp3 = bsl.t10;


    tmp7 = _mm_slli_epi32(tmp3, 31);
    tmp8 = _mm_slli_epi32(tmp3, 30);
    tmp9 = _mm_slli_epi32(tmp3, 25);
    tmp7 = _mm_xor_si128(tmp7, tmp8);
    tmp7 = _mm_xor_si128(tmp7, tmp9);
    tmp8 = _mm_srli_si128(tmp7, 4);
    tmp7 = _mm_slli_si128(tmp7, 12);
    tmp3 = _mm_xor_si128(tmp3, tmp7);
    tmp2 = _mm_srli_epi32(tmp3, 1);
    tmp4 = _mm_srli_epi32(tmp3, 2);
    tmp5 = _mm_srli_epi32(tmp3, 7);
    tmp2 = _mm_xor_si128(tmp2, tmp4);
    tmp2 = _mm_xor_si128(tmp2, tmp5);
    tmp2 = _mm_xor_si128(tmp2, tmp8);
    tmp3 = _mm_xor_si128(tmp3, tmp2);
    tmp6 = _mm_xor_si128(tmp6, tmp3);

    return tmp6;
}

__m128i gfmul_original_docA_raw(__m128i a, __m128i b){
    __m128i tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9;
    tmp3 = _mm_clmulepi64_si128(a, b, 0x00);
    tmp4 = _mm_clmulepi64_si128(a, b, 0x10);
    tmp5 = _mm_clmulepi64_si128(a, b, 0x01);
    tmp6 = _mm_clmulepi64_si128(a, b, 0x11);
    tmp4 = _mm_xor_si128(tmp4, tmp5);
    tmp5 = _mm_slli_si128(tmp4, 8);
    tmp4 = _mm_srli_si128(tmp4, 8);
    tmp3 = _mm_xor_si128(tmp3, tmp5);
    tmp6 = _mm_xor_si128(tmp6, tmp4);
    tmp7 = _mm_srli_epi32(tmp3, 31);
    tmp8 = _mm_srli_epi32(tmp6, 31);
    tmp3 = _mm_slli_epi32(tmp3, 1);
    tmp6 = _mm_slli_epi32(tmp6, 1);
    tmp9 = _mm_srli_si128(tmp7, 12);
    tmp8 = _mm_slli_si128(tmp8, 4);
    tmp7 = _mm_slli_si128(tmp7, 4);
    tmp3 = _mm_or_si128(tmp3, tmp7);
    tmp6 = _mm_or_si128(tmp6, tmp8);
    tmp6 = _mm_or_si128(tmp6, tmp9);
    tmp7 = _mm_slli_epi32(tmp3, 31);
    tmp8 = _mm_slli_epi32(tmp3, 30);
    tmp9 = _mm_slli_epi32(tmp3, 25);
    tmp7 = _mm_xor_si128(tmp7, tmp8);
    tmp7 = _mm_xor_si128(tmp7, tmp9);
    tmp8 = _mm_srli_si128(tmp7, 4);
    tmp7 = _mm_slli_si128(tmp7, 12);
    tmp3 = _mm_xor_si128(tmp3, tmp7);
    tmp2 = _mm_srli_epi32(tmp3, 1);
    tmp4 = _mm_srli_epi32(tmp3, 2);
    tmp5 = _mm_srli_epi32(tmp3, 7);
    tmp2 = _mm_xor_si128(tmp2, tmp4);
    tmp2 = _mm_xor_si128(tmp2, tmp5);
    tmp2 = _mm_xor_si128(tmp2, tmp8);
    tmp3 = _mm_xor_si128(tmp3, tmp2);
    tmp6 = _mm_xor_si128(tmp6, tmp3);

    return tmp6;
}


/*
 * Testing code:
 *  __m128i x1 = _mm_set_epi32(0x00000000, 0x00000010, 0x00000000, 0x00000010);
    __m128i x2 = _mm_set_epi32(0x00000000, 0x00000020, 0x00000000, 0x00000020);
    struct int256 xRet = bitshift_left256(x2, x1, 3);
    qInfo() <<"[x2,x1] << 3: "<<print128_hex_lanes(xRet.t32)<<" "<<print128_hex_lanes(xRet.t10);
    xRet = bitshift_left256(x2, x1, 65);
    qInfo() <<"[x2,x1] << 65: "<<print128_hex_lanes(xRet.t32)<<" "<<print128_hex_lanes(xRet.t10);
 *
 */
struct int256 bitshift_left256(__m128i i32, __m128i i10, unsigned char count){
    struct int256 ret;
    if (count >= 128){
        ret.t10 = _mm_setzero_si128();
        ret.t32 = bitshift_left(i10, count-128);
        return ret;
    }
    // else
    __m128i carry = bitshift_right(i10, 128-count);
    ret.t10 = bitshift_left(i10, count);
    ret.t32 = _mm_or_si128(bitshift_left(i32, count), carry);

    return ret;
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

__m128i bitshift_right(__m128i a, unsigned char count){
    __m128i carry = _mm_srli_si128(a, 8);
    if (count >= 64){
        return _mm_srli_epi64(carry, count-64);
    }
    //else
    carry = _mm_slli_epi64(carry, 64-count);

    a = _mm_srli_epi64(a, count);
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

QString print256_hex_lanes(struct int256 var)
{
    uint32_t valHigh[4];
    memcpy(valHigh, &(var.t32), sizeof(valHigh));
    uint32_t valLow[4];
    memcpy(valLow, &(var.t10), sizeof(valLow));

    // print high lane first
    QString s = QString::asprintf(
        "%08X %08X %08X %08X %08X %08X %08X %08X\n",
        valHigh[3], valHigh[2], valHigh[1], valHigh[0], valLow[3], valLow[2], valLow[1], valLow[0]);
    return s;
}

#include <QElapsedTimer>
#include "gcm/benchmarkutil.h"
void gfmul_test(){
    qInfo() << "hasSSSE3?: " << hasSSSE3();
    /* Test vektoren (TEST 1 - nach Intel Doc 2014 - Page 78 in Doc A.: "Intel Carry-Less Multiplication Instruction and its Usage for Computing in GCM mode"):
    / a: 7b5b5465 73745665 63746f72 5d53475d
    / b: 48692853 68617929 5b477565 726f6e5d
    / q: 00000000 00000000 00000000 00000087
    / assert: GFMUL(a,b,q) = c = 040229a0 9a5ed12e 7e4e10da 323506d2
    */
    __m128i a = _mm_set_epi32(0x7b5b5465, 0x73745665, 0x63746f72, 0x5d53475d);
    __m128i b = _mm_set_epi32(0x48692853, 0x68617929, 0x5b477565, 0x726f6e5d);
    __m128i res_assert = _mm_set_epi32(0x040229a0, 0x9a5ed12e, 0x7e4e10da, 0x323506d2);
    __m128i q = _mm_set_epi32(0x00000000, 0x00000000, 0x00000000, 0x00000087);
    qInfo() << "a: "<<print128_hex_lanes(a)<<", b: "<<print128_hex_lanes(b);
    __m128i res = gfmul_k_optimized(a,b,q);
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
    q = _mm_set_epi32(0x00000000, 0x00000000, 0x00000000, 0x00000087);
    __m128i a_refl = reflect_xmm(a);
    __m128i b_refl = reflect_xmm(b);
    __m128i c_assert = _mm_set_epi32(0x065B7FC3, 0x340123F2, 0x6DDAA34B, 0x50D7CA5B);
    __m128i c_refl_assert = _mm_set_epi32(0xda53eb0a, 0xd2c55bb6, 0x4fc4802c, 0xc3feda60);
    qInfo() << "a: "<<print128_hex_lanes(a)<<", b: "<<print128_hex_lanes(b);
    qInfo() << "a_refl: "<<print128_hex_lanes(a_refl)<<", b_refl: "<<print128_hex_lanes(b_refl);
    //qInfo() << "q_refl: "<<print128_hex_lanes(q_refl);
    __m128i c = gfmul_k_optimized(a_refl, b_refl, q);       // <<< switch here between gfmul and gfmul_k_optimized
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
    res = gfmul_reversed_bl_opt(a,b);
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


    /**
     * Timing test to identify whether DocA or DocB algo is better
     *
     */
    BenchmarkUtil::run("DocA-raw", gfmul_original_docA_raw, a, b);
    BenchmarkUtil::run("DocA-mod", gfmul_original_docA_mod, a, b);
    BenchmarkUtil::run("DocB", gfmul_reversed, a, b);
    BenchmarkUtil::run("DocB-BlOpt", gfmul_reversed_bl_opt, a, b);

    /*
     * Test of gfmul_reversed_k_optimized
     *
     * We try to get around using reflect_xmm, thats why OUT should already be the final value
     *
     * Assert against res: da53eb0a d2c55bb6 4fc4802c c3feda60
     * res_refl is: 065B7FC3 340123F2 6DDAA34B 50D7CA5B
     *
     */
    qInfo() <<"\nNow testing gfmul_reversed_k_optimized";
    a = _mm_set_epi32(0x952b2a56, 0xa5604ac0, 0xb32b6656, 0xa05b40b6);
    b = _mm_set_epi32(0xdfa6bf4d, 0xed81db03, 0xffcaff95, 0xf830f061);
    res_assert = _mm_set_epi32(0xda53eb0a, 0xd2c55bb6, 0x4fc4802c, 0xc3feda60);
    qInfo() << "a: "<<print128_hex_lanes(a)<<", b: "<<print128_hex_lanes(b);
    res = gfmul_reversed_k_optimized(a,b);
    qInfo() << "res: (a, b, q):\n|>"<<print128_hex_lanes(res);
    qInfo() << "res_assert:\n|>"<<print128_hex_lanes(res_assert);
    if(_mm_test_all_zeros(_mm_set_epi32(0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff), _mm_xor_si128(res_assert, res))) {
        qInfo() << "Assertion (res == res_assert): holds true!";
    } else {
        qWarning() << "Assertion (res == res_assert): is false!";
    }
    qInfo("\n");


    /*
     * Test of gfmul_k_optimized_reversed
     *
     * This is computing a reversed outcome, e.g. c from TEST-2
     * Assert against c: 065B7FC3 340123F2 6DDAA34B 50D7CA5B
     *
     */
    qInfo() <<"\nNow testing gfmul_k_optimized_reversed";
    a = _mm_set_epi32(0x952b2a56, 0xa5604ac0, 0xb32b6656, 0xa05b40b6);
    b = _mm_set_epi32(0xdfa6bf4d, 0xed81db03, 0xffcaff95, 0xf830f061);
    q = _mm_set_epi32(0x00000000, 0x00000000, 0x00000000, 0x00000087);
    res_assert = _mm_set_epi32(0x065B7FC3, 0x340123F2, 0x6DDAA34B, 0x50D7CA5B);
    qInfo() << "a: "<<print128_hex_lanes(a)<<", b: "<<print128_hex_lanes(b);
    res = gfmul_k_optimized_reversed(a,b, q);
    res_refl = reflect_xmm(res);
    qInfo() << "res: (a, b, q):\n|>"<<print128_hex_lanes(res);
    qInfo() << "res_refl: \n|>"<<print128_hex_lanes(res_refl);
    if(_mm_test_all_zeros(_mm_set_epi32(0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff), _mm_xor_si128(res_assert, res))) {
        qInfo() << "Assertion (res == res_assert): holds true!";
    } else {
        qWarning() << "Assertion (res == res_assert): is false!";
    }
    qInfo("\n");
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
