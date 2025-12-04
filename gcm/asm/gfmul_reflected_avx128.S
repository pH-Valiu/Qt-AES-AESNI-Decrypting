#--------------------
# CONSTANTS
#--------------------

.section .rodata

# 64 byte alignment (2^4)
#ifdef __ELF__
.p2align 4
#else
.balign 16
#endif

.global Q_r_128_asm

#ifdef __ELF__
.type Q_r_128_asm, @object
#endif

# this is equivalent to _mm_set_epi32(0, 0, 0xc2000000, 0);
Q_r_128_asm:
    .quad 0xc200000000000000
    .quad 0x0000000000000000

#ifdef __ELF__
.size Q_r_128_asm, 16   # tell the linker its size
#endif


#--------------------
# ALGORITHM
#--------------------

.text
.intel_syntax noprefix
.global gfmul_reflected_avx128

#ifdef __ELF__
.type gfmul_reflected_avx128, @function
#endif

gfmul_reflected_avx128:
    # Parameters are passed through pointers: RCX -> a, RDX -> b, R8 -> result:
    # xmm0 := a
    # xmm1 := b

    # since all xmm registers are caller saved, we can use all 16 (32 with AVX-512 support) registers
    # I want to avoid using the stack
    movdqu xmm0, [rcx]
    movdqu xmm1, [rdx]

    # registers xmm2 to xmm5 store all xmm0 = a
    movdqa xmm2, xmm0
    movdqa xmm3, xmm0
    movdqa xmm4, xmm0
    movdqa xmm5, xmm0

    pclmulqdq xmm2, xmm1, 0x00   # xmm2 := a0b0
    pclmulqdq xmm3, xmm1, 0x01   # xmm3 := a1b0
    pclmulqdq xmm4, xmm1, 0x10   # xmm4 := a0b1
    pclmulqdq xmm5, xmm1, 0x11   # xmm5 := a1b1

    pxor xmm3, xmm4             # xmm3 := mid = a1b0 + a0b1
    movdqa xmm7, xmm3           # xmm7 := xmm3 = copy(mid)

    pslldq xmm3, 8              # xmm3 := mid <<_epi128 64
    pxor xmm2, xmm3             # xmm2 := xmm2 xor xmm3 = lo =  a0b0 xor (mid <<_epi128 64)

    psrldq xmm7, 8              # xmm7 := mid >>_epi128 64
    pxor xmm5, xmm7             # xmm5 := xmm5 xor xmm7 = hi = a1b1 xor (mid >>_epi128 64)

    # Now xmm5 := hi; xmm2 := lo
    movdqa xmm0, xmm5           # xmm0 := xmm5 = copy(hi)
    movdqa xmm1, xmm2           # xmm1 := xmm2 = copy(lo)

    # Step 1.1: Bitshift <<_epi256 1
    psrld xmm1, 31               # xmm1 := lo_32bitMSB = lo >>_epi32 31
    psrld xmm0, 31               # xmm0 := hi_32bitMSB = hi >>_epi32 31
    pslld xmm2, 1                # xmm2 := lo = lo <<_epi32 1
    pslld xmm5, 1                # xmm5 := hi = hi <<_epi32 1
    movdqa xmm6, xmm1            # xmm6 := xmm1 = copy(lo_32bitMSB)
    psrldq xmm6, 12              # xmm6 := lo_128bitMSB = lo_32bitMSB >>_epi128 96
    pslldq xmm1, 4               # xmm1 := lo_32bitMSB = lo_32bitMSB <<_epi128 32
    pslldq xmm0, 4               # xmm0 := hi_32bitMSB = hi_32bitMSB <<_epi128 32
    por xmm2, xmm1               # xmm2 := lo = lo | lo_32bitMSB
    por xmm5, xmm0               # xmm5 := hi = hi | hi_32bitMSB
    por xmm5, xmm6               # xmm5 := hi = hi | lo_128bitMSB

    # Now xmm5 := hi; xmm2 := lo
    movdqa xmm0, xmm2           # xmm0 := xmm2 = copy(lo)
    movdqa xmm1, xmm2           # xmm1 := xmm2 = copy(lo)

    # Step 2.0: Load constants
    movdqu xmm6, [rip + Q_r_128_asm]     # xmm6 := Q_r_128
    pxor xmm8, xmm8                      # xmm8 := ZERO (storing ZERO as constant and loading twice is more load than just doing one xor every time)

    # Step 2.1: Reduce
    pclmulqdq xmm1, xmm6, 0x00          # xmm1 := x = lo[0] * Q_r[0]
    movdqa xmm7, xmm1                   # xmm7 := xmm1 = copy(x)
    psrldq xmm7, 8                      # xmm7 := x >>_epi128 64
    pxor xmm5, xmm7                     # xmm5 := hi = hi xor (x >>_epi128 64)
    punpcklqdq xmm0, xmm8               # xmm0 := _mm_unpacklo_epi64(lo, ZERO)
    pxor xmm5, xmm0                     # xmm5 := hi = hi xor (_mm512_unpacklo_epi64(lo, ZERO))
    pslldq xmm1, 8                      # xmm1 := x <<_epi128 64
    pxor xmm2, xmm1                     # xmm2 := lo = lo xor (x <<_epi128 64)

    # Step 2.2: Reduce
    pclmulqdq xmm6, xmm2, 0x10          # xmm6 := x = lo[1] * Q_r[0]
    pxor xmm5, xmm6                     # xmm5 := hi = hi xor x
    punpckhqdq xmm8, xmm2               # xmm8 := _mm_unpackhi_epi64(ZERO, lo)
    pxor xmm5, xmm8                     # xmm5 := hi = hi xor (_mm512_unpackhi_epi64(ZERO, lo))


    # Step 3: Return hi
    # Return value must be supplied through:
    # zmm0 := returnVal
    # or into pointer
    movdqu [r8], xmm5
    ret

#ifdef __ELF__
.size gfmul_reflected_avx128, .-gfmul_reflected_avx128
#endif

.att_syntax
