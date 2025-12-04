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
    .quad 0x0000000000000000
    .quad 0xc200000000000000
    .quad 0x0000000000000000
    .quad 0xc200000000000000
    .quad 0x0000000000000000
    .quad 0xc200000000000000
    .quad 0x0000000000000000
    .quad 0xc200000000000000

#ifdef __ELF__
.size Q_r_512_asm, 64   # tell the linker its size
#endif


#--------------------
# ALGORITHM
#--------------------

.text
.intel_syntax noprefix
.global gfmul_reflected_avx512

#ifdef __ELF__
.type gfmul_reflected_avx512, @function
#endif

gfmul_reflected_avx512:
    # Parameters are passed through registers:
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

#ifdef __ELF__
.size gfmul_reflected_avx512, .-gfmul_reflected_avx512
#endif

.att_syntax
