#--------------------
# CONSTANTS
#--------------------

.section .rodata

# 64 byte alignment (2^6)
#ifdef __ELF__
.p2align 6
#else
.balign 64
#endif

.global Q_r_512_asm

#ifdef __ELF__
.type Q_r_512_asm, @object
#endif

Q_r_512_asm:
    .quad 0xc200000000000000
    .quad 0x0000000000000000
    .quad 0xc200000000000000
    .quad 0x0000000000000000
    .quad 0xc200000000000000
    .quad 0x0000000000000000
    .quad 0xc200000000000000
    .quad 0x0000000000000000

#ifdef __ELF__
.size Q_r_512_asm, 64   # tell the linker its size
#endif


#--------------------
# ALGORITHM
#--------------------

.text
.intel_syntax noprefix
.global avx512_xor_si512_asm
.global gfmul_reflected_avx512_512
.global gfmul_reflected_avx512_128
.global gfmul_reflected_avx512_512_test

#ifdef __ELF__
.type gfmul_reflected_avx512_512, @function
.type gfmul_reflected_avx512_128, @function
.type gfmul_reflected_avx512_512_test, @function
#endif

avx512_xor_si512_asm:
    # expect rcx -> &a, rdx -> &b
    # expect rcx -> &result
    vmovdqu64 zmm0, [rcx]
    vmovdqu64 zmm1, [rdx]

    vpxord zmm0, zmm0, zmm1

    vmovdqu64 [rcx], zmm0
    ret
# This is faster than gfmul_reflected_avx128 because it uses AVX-512 instructions
gfmul_reflected_avx512_128:
    vmovdqu xmm0, [rcx]
    vmovdqu xmm1, [rdx]

    vpclmulqdq xmm2, xmm0, xmm1, 0x00   # a0b0
    vpclmulqdq xmm3, xmm0, xmm1, 0x01   # a1b0
    vpclmulqdq xmm4, xmm0, xmm1, 0x10   # a0b1
    vpclmulqdq xmm5, xmm0, xmm1, 0x11   # a1b1

    vpxord xmm3, xmm3, xmm4             # xmm3 := mid = a1b0 + a0b1

    vpslldq xmm4, xmm3, 8               # xmm4 := mid <<_epi128 64
    vpxord xmm1, xmm2, xmm4             # xmm1 := lo = a0b0 xor (mid <<_epi128 64)

    vpsrldq xmm4, xmm3, 8               # xmm4 := mid >>_epi128 64
    vpxord xmm0, xmm5, xmm4             # xmm0 := hi = a1b1 xor (mid >>_epi128 64)

    # Now xmm0 := hi; xmm1 := lo

    # Step 1.1: Bitshift <<_epi1024 1
    vpsrld xmm2, xmm1, 31               # xmm2 := lo_32bitMSB = lo >>_epi32 31
    vpsrld xmm3, xmm0, 31               # xmm3 := hi_32bitMSB = hi >>_epi32 31
    vpslld xmm1, xmm1, 1                # xmm1 := lo = lo <<_epi32 1
    vpslld xmm0, xmm0, 1                # xmm0 := hi = hi <<_epi32 1
    vpsrldq xmm4, xmm2, 12              # xmm4 := lo_128bitMSB = lo_32bitMSB >>_epi128 96
    vpslldq xmm2, xmm2, 4               # xmm2 := lo_32bitMSB = lo_32bitMSB <<_epi128 32
    vpslldq xmm3, xmm3, 4               # xmm3 := hi_32bitMSB = hi_32bitMSB <<_epi128 32
    vpord xmm1, xmm1, xmm2              # xmm1 := lo = lo | lo_32bitMSB
    vpord xmm0, xmm0, xmm3              # xmm0 := hi = hi | hi_32bitMSB
    vpord xmm0, xmm0, xmm4              # xmm0 := hi = hi | lo_128bitMSB

    # Step 2.0: Load constants
    vmovdqu xmm6, [rip + Q_r_128_asm]     # xmm6 := Q_r_128
    vpxord xmm7, xmm7, xmm7                 # xmm7 := ZERO (storing ZERO as constant and loading twice is more load than just doing one xor every time)

    # Step 2.1: Reduce
    vpclmulqdq xmm2, xmm1, xmm6, 0x00   # xmm2 := x = lo_epi128[0] * Q_r_epi128[0]
    vpsrldq xmm3, xmm2, 8               # xmm3 := x >>_epi128 64
    vpxord xmm0, xmm0, xmm3             # xmm0 := hi = hi xor (x >>_epi128 64)
    vpunpcklqdq xmm3, xmm1, xmm7        # xmm3 := _mm_unpacklo_epi64(lo, ZERO)
    vpxord xmm0, xmm0, xmm3             # xmm0 := hi = hi xor (_mm_unpacklo_epi64(lo, ZERO))
    vpslldq xmm3, xmm2, 8               # xmm3 := x <<_epi128 64
    vpxord xmm1, xmm1, xmm3             # xmm1 := lo = lo xor (x <<_epi128 64)

    # Step 2.2: Reduce
    vpclmulqdq xmm2, xmm1, xmm6, 0x01   # xmm2 := x = lo_epi128[1] * Q_r_epi128[0]
    vpxord xmm0, xmm0, xmm2             # xmm0 := hi = hi xor x
    vpunpckhqdq xmm3, xmm7, xmm1        # xmm3 := _mm_unpackhi_epi64(ZERO, lo)
    vpxord xmm0, xmm0, xmm3             # xmm0 := hi = hi xor (_mm_unpackhi_epi64(ZERO, lo))


    # Step 3: Return hi
    # Return value must be supplied through:
    # xmm0 := returnVal
    vmovdqu [r8], xmm0
    ret



# The goal is to perform gfmul_reflected on all 4 quad words (__m128i) at the same time
# Therefore, each step must be seen as happening independept on each quad word
# To understand the commands, look at the C++ implementation of gfmul_reflected for __m128i and transpose to __m512

# The main reason for implementing it in assembly directly is because of a known old MinGW Windows bug that causes alignment issues when off shoring __m512 to stack.
# Therefore, __m512 representation can only be used temporary for specific computations.
# For the rest, we must use assembler code and pass parameters and return buffer as pointers.

# registers used:
# zmm0,zmm1,zmm2,zmm3,zmm4,zmm5,zmm6,zmm7
gfmul_reflected_avx512_512:
    # Parameters are passed through pointers: RCX -> a, RDX -> b, R8 -> result:
    # zmm0 := a
    # zmm1 := b

    # since all zmm registers are caller saved, we can use all 32 registers
    # I want to avoid using the stack

    # Step 1.0: Multiplication
    vmovdqu64 zmm0, [rcx]
    vmovdqu64 zmm1, [rdx]

    vpclmulqdq zmm2, zmm0, zmm1, 0x00   # a0b0
    vpclmulqdq zmm3, zmm0, zmm1, 0x01   # a1b0
    vpclmulqdq zmm4, zmm0, zmm1, 0x10   # a0b1
    vpclmulqdq zmm5, zmm0, zmm1, 0x11   # a1b1

    vpxord zmm3, zmm3, zmm4             # zmm3 := mid = a1b0 + a0b1

    vpslldq zmm4, zmm3, 8               # zmm4 := mid <<_epi128 64
    vpxord zmm1, zmm2, zmm4             # zmm1 := lo = a0b0 xor (mid <<_epi128 64)

    vpsrldq zmm4, zmm3, 8               # zmm4 := mid >>_epi128 64
    vpxord zmm0, zmm5, zmm4             # zmm0 := hi = a1b1 xor (mid >>_epi128 64)

    # Now zmm0 := hi; zmm1 := lo

    # Step 1.1: Bitshift <<_epi1024 1
    vpsrld zmm2, zmm1, 31               # zmm2 := lo_32bitMSB = lo >>_epi32 31
    vpsrld zmm3, zmm0, 31               # zmm3 := hi_32bitMSB = hi >>_epi32 31
    vpslld zmm1, zmm1, 1                # zmm1 := lo = lo <<_epi32 1
    vpslld zmm0, zmm0, 1                # zmm0 := hi = hi <<_epi32 1
    vpsrldq zmm4, zmm2, 12              # zmm4 := lo_128bitMSB = lo_32bitMSB >>_epi128 96
    vpslldq zmm2, zmm2, 4               # zmm2 := lo_32bitMSB = lo_32bitMSB <<_epi128 32
    vpslldq zmm3, zmm3, 4               # zmm3 := hi_32bitMSB = hi_32bitMSB <<_epi128 32
    vpord zmm1, zmm1, zmm2              # zmm1 := lo = lo | lo_32bitMSB
    vpord zmm0, zmm0, zmm3              # zmm0 := hi = hi | hi_32bitMSB
    vpord zmm0, zmm0, zmm4              # zmm0 := hi = hi | lo_128bitMSB

    # Step 2.0: Load constants
    vmovdqu64 zmm6, [rip + Q_r_512_asm]     # zmm6 := Q_r_512
    vpxord zmm7, zmm7, zmm7                 # zmm7 := ZERO (storing ZERO as constant and loading twice is more load than just doing one xor every time)

    # Step 2.1: Reduce
    vpclmulqdq zmm2, zmm1, zmm6, 0x00   # zmm2 := x = lo_epi128[0] * Q_r_epi128[0]
    vpsrldq zmm3, zmm2, 8               # zmm3 := x >>_epi128 64
    vpxord zmm0, zmm0, zmm3             # zmm0 := hi = hi xor (x >>_epi128 64)
    vpunpcklqdq zmm3, zmm1, zmm7        # zmm3 := _mm512_unpacklo_epi64(lo, ZERO)
    vpxord zmm0, zmm0, zmm3             # zmm0 := hi = hi xor (_mm512_unpacklo_epi64(lo, ZERO))
    vpslldq zmm3, zmm2, 8               # zmm3 := x <<_epi128 64
    vpxord zmm1, zmm1, zmm3             # zmm1 := lo = lo xor (x <<_epi128 64)

    # Step 2.2: Reduce
    vpclmulqdq zmm2, zmm1, zmm6, 0x01   # zmm2 := x = lo_epi128[1] * Q_r_epi128[0]
    vpxord zmm0, zmm0, zmm2             # zmm0 := hi = hi xor x
    vpunpckhqdq zmm3, zmm7, zmm1        # zmm3 := _mm512_unpackhi_epi64(ZERO, lo)
    vpxord zmm0, zmm0, zmm3             # zmm0 := hi = hi xor (_mm512_unpackhi_epi64(ZERO, lo))


    # Step 3: Return hi
    # Return value must be supplied through:
    # zmm0 := returnVal
    vmovdqu64 [r8], zmm0
    ret

gfmul_reflected_avx512_512_test:
    # Parameters are passed through pointers: RCX -> a, RDX -> b, R8 -> result:
    # zmm0 := a
    # zmm1 := b

    # since all zmm registers are caller saved, we can use all 32 registers
    # I want to avoid using the stack

    # Step 1.0: Multiplication
    vmovdqu xmm0, [rcx]
    vmovdqu xmm1, [rdx]
    vmovdqu xmm2, [r8]
    vmovdqu xmm3, [r9]
    mov rcx, [rsp + 40]
    mov rdx, [rsp + 48]
    mov r8, [rsp + 56]
    mov r9, [rsp + 64]
    vmovdqu xmm4, [rcx]
    vmovdqu xmm5, [rdx]
    vmovdqu xmm6, [r8]
    vmovdqu xmm7, [r9]

    vinserti64x2 zmm0, zmm0, xmm0, 0
    vinserti64x2 zmm0, zmm0, xmm1, 1
    vinserti64x2 zmm0, zmm0, xmm2, 2
    vinserti64x2 zmm0, zmm0, xmm3, 3

    vinserti64x2 zmm1, zmm1, xmm4, 0
    vinserti64x2 zmm1, zmm1, xmm5, 1
    vinserti64x2 zmm1, zmm1, xmm6, 2
    vinserti64x2 zmm1, zmm1, xmm7, 3

    vpclmulqdq zmm2, zmm0, zmm1, 0x00   # a0b0
    vpclmulqdq zmm3, zmm0, zmm1, 0x01   # a1b0
    vpclmulqdq zmm4, zmm0, zmm1, 0x10   # a0b1
    vpclmulqdq zmm5, zmm0, zmm1, 0x11   # a1b1

    vpxord zmm3, zmm3, zmm4             # zmm3 := mid = a1b0 + a0b1

    vpslldq zmm4, zmm3, 8               # zmm4 := mid <<_epi128 64
    vpxord zmm1, zmm2, zmm4             # zmm1 := lo = a0b0 xor (mid <<_epi128 64)

    vpsrldq zmm4, zmm3, 8               # zmm4 := mid >>_epi128 64
    vpxord zmm0, zmm5, zmm4             # zmm0 := hi = a1b1 xor (mid >>_epi128 64)

    # Now zmm0 := hi; zmm1 := lo

    # Step 1.1: Bitshift <<_epi1024 1
    vpsrld zmm2, zmm1, 31               # zmm2 := lo_32bitMSB = lo >>_epi32 31
    vpsrld zmm3, zmm0, 31               # zmm3 := hi_32bitMSB = hi >>_epi32 31
    vpslld zmm1, zmm1, 1                # zmm1 := lo = lo <<_epi32 1
    vpslld zmm0, zmm0, 1                # zmm0 := hi = hi <<_epi32 1
    vpsrldq zmm4, zmm2, 12              # zmm4 := lo_128bitMSB = lo_32bitMSB >>_epi128 96
    vpslldq zmm2, zmm2, 4               # zmm2 := lo_32bitMSB = lo_32bitMSB <<_epi128 32
    vpslldq zmm3, zmm3, 4               # zmm3 := hi_32bitMSB = hi_32bitMSB <<_epi128 32
    vpord zmm1, zmm1, zmm2              # zmm1 := lo = lo | lo_32bitMSB
    vpord zmm0, zmm0, zmm3              # zmm0 := hi = hi | hi_32bitMSB
    vpord zmm0, zmm0, zmm4              # zmm0 := hi = hi | lo_128bitMSB

    # Step 2.0: Load constants
    vmovdqu64 zmm6, [rip + Q_r_512_asm]     # zmm6 := Q_r_512
    vpxord zmm7, zmm7, zmm7                 # zmm7 := ZERO (storing ZERO as constant and loading twice is more load than just doing one xor every time)

    # Step 2.1: Reduce
    vpclmulqdq zmm2, zmm1, zmm6, 0x00   # zmm2 := x = lo_epi128[0] * Q_r_epi128[0]
    vpsrldq zmm3, zmm2, 8               # zmm3 := x >>_epi128 64
    vpxord zmm0, zmm0, zmm3             # zmm0 := hi = hi xor (x >>_epi128 64)
    vpunpcklqdq zmm3, zmm1, zmm7        # zmm3 := _mm512_unpacklo_epi64(lo, ZERO)
    vpxord zmm0, zmm0, zmm3             # zmm0 := hi = hi xor (_mm512_unpacklo_epi64(lo, ZERO))
    vpslldq zmm3, zmm2, 8               # zmm3 := x <<_epi128 64
    vpxord zmm1, zmm1, zmm3             # zmm1 := lo = lo xor (x <<_epi128 64)

    # Step 2.2: Reduce
    vpclmulqdq zmm2, zmm1, zmm6, 0x01   # zmm2 := x = lo_epi128[1] * Q_r_epi128[0]
    vpxord zmm0, zmm0, zmm2             # zmm0 := hi = hi xor x
    vpunpckhqdq zmm3, zmm7, zmm1        # zmm3 := _mm512_unpackhi_epi64(ZERO, lo)
    vpxord zmm0, zmm0, zmm3             # zmm0 := hi = hi xor (_mm512_unpackhi_epi64(ZERO, lo))


    # Step 3: Return hi
    # Return value must be supplied through:
    # zmm0 := returnVal
    mov rcx, [rsp + 72]
    mov rdx, [rsp + 80]
    mov r8, [rsp + 88]
    mov r9, [rsp + 96]

    vextracti64x2 xmm3, zmm0, 0
    vextracti64x2 xmm4, zmm0, 1
    vextracti64x2 xmm5, zmm0, 2
    vextracti64x2 xmm6, zmm0, 3

    vmovdqu [rcx], xmm3
    vmovdqu [rdx], xmm4
    vmovdqu [r8], xmm5
    vmovdqu [r9], xmm6
    ret

    #vmovdqu64 [r8], zmm0
    #ret


#ifdef __ELF__
.size gfmul_reflected_avx512_512, .-gfmul_reflected_avx512_512
.size gfmul_reflected_avx512_128, .-gfmul_reflected_avx512_128
.size gfmul_reflected_avx512_512_test, .-gfmul_reflected_avx512_512_test
#endif

.att_syntax
