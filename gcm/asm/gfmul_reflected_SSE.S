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
# This algorithm is suited for an environment with just SSE instructions without AVX-512

.text
.intel_syntax noprefix
.global gfmul_reflected_avx128
.global gfmul_reflected_avx128_times_four

#ifdef __ELF__
.type gfmul_reflected_avx128, @function
.type gfmul_reflected_avx128_times_four, @function
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
    # xmm0 := returnVal
    # or into pointer
    movdqu [r8], xmm5
    ret

gfmul_reflected_avx128_times_four:
    # this method takes 9 parameters (all pointers to __m128i values)
    # param 1-4: Y1, Y2, Y3, Y4
    # param 5-8: H1, H2, H3, H4
    # param   9: out
    #
    # because of Windows x86 System V ABI, only the first four arguments are passed through registers [rcx, rdx, r8, r9]
    # the remaining ones are pushed onto the stack such that it looks the following way
    # --- stack HIGH ADDRESSES ---
    # arg &out
    # arg &H4
    # arg &H3
    # arg &H2
    # arg &H1
    # [shadow space for &Y4 - left empty]
    # [shadow space for &Y3 - left empty]
    # [shadow space for &Y2 - left empty]
    # [shadow space for &Y1 - left empty]
    # returnAddress     << rsp is looking here down
    # --- stack LOW ADDRESSES ---

    # -- load parameters into xmm0 - xmm7 --
    movdqu xmm0, [rcx]  # xmm0 := Y1
    movdqu xmm1, [rdx]  # xmm1 := Y2
    movdqu xmm2, [r8]   # xmm2 := Y3
    movdqu xmm3, [r9]   # xmm3 := Y4

    mov rcx, [rsp + 40] # rcx := &H1
    mov rdx, [rsp + 48] # rdx := &H2
    mov r8, [rsp + 56]  # r8 := &H3
    mov r9, [rsp + 64]  # r9 := &H4
    movdqu xmm4, [rcx]  # xmm4 := H1
    movdqu xmm5, [rdx]  # xmm5 := H2
    movdqu xmm6, [r8]   # xmm6 := H3
    movdqu xmm7, [r9]   # xmm7 := H4

    # --- Step 1: Karatsuba multiplication to yield [C3:C2:C1:C0]

    # make copys of H1-H4
    movdqu xmm8, xmm4   # xmm8 := H1 copy
    movdqu xmm9, xmm5   # xmm9 := H2 copy
    movdqu xmm10, xmm6  # xmm10 := H3 copy
    movdqu xmm11, xmm7  # xmm11 := H4 copy

    # -- Step 1.1: compute C0 of [C3:C2:C1:C0] --
    pclmulqdq xmm8, xmm0, 0x00  # xmm8 := H1_Y1_lo
    pclmulqdq xmm9, xmm1, 0x00  # xmm9 := H2_Y2_lo
    pclmulqdq xmm10, xmm2, 0x00 # xmm10 := H3_Y3_lo
    pclmulqdq xmm11, xmm3, 0x00 # xmm11 := H4_Y4_lo
    pxor xmm8, xmm9             # xmm8 := xmm8 xor xmm9
    pxor xmm8, xmm10            # xmm8 := xmm8 xor xmm10
    pxor xmm8, xmm11            # xmm8 := xmm8 xor xmm11

    # make copys of H1-H4 again
    movdqu xmm9, xmm4           # xmm9 := H1 copy
    movdqu xmm10, xmm5          # xmm10 := H2 copy
    movdqu xmm11, xmm6          # xmm11 := H3 copy
    movdqu xmm12, xmm7          # xmm12 := H4 copy

    # -- Step 1.2: compute C3 of [C3:C2:C1:C0]
    pclmulqdq xmm9, xmm0, 0x11  # xmm9 := H1_Y1_hi
    pclmulqdq xmm10, xmm1, 0x11 # xmm10 := H2_Y2_hi
    pclmulqdq xmm11, xmm2, 0x11 # xmm11 := H3_Y3_hi
    pclmulqdq xmm12, xmm3, 0x11 # xmm12 := H4_Y4_hi
    pxor xmm9, xmm10            # xmm9 := xmm9 xor xmm10
    pxor xmm9, xmm11            # xmm9 := xmm9 xor xmm11
    pxor xmm9, xmm12            # xmm9 := xmm9 xor xmm12


    # --> xmm8 := C0, xmm9 := C3

    # -- Step 1.3: compute (H1_high xor H1_low) CLMUL (Y1_high xor Y1_low) for each pair respectively

    # - xmm0 := H1_Y1_mid
    pshufd xmm10, xmm4, 78  # xmm10 := _mm_shuffle_epi32(H1, 78)
    pshufd xmm11, xmm0, 78  # xmm11 := _mm_shuffle_epi32(Y1, 78)
    pxor xmm4, xmm10        # xmm4 := xmm4 xor xmm10 = H1 xor _mm_shuffle_epi32(H1, 78)
    pxor xmm0, xmm11        # xmm0 := xmm0 xor xmm11 = Y1 xor _mm_shuffle_epi32(Y1, 78)
    pclmulqdq xmm0, xmm4, 0x00  # xmm0 := middle term of H1_Y1

    # - xmm1 := H2_Y2_mid
    pshufd xmm10, xmm5, 78      # xmm10 := shuffle(H2, 78)
    pshufd xmm11, xmm1, 78      # xmm11 := shuffle(Y2, 78)
    pxor xmm5, xmm10            # xmm5 := H2 xor shuffle(H2, 78)
    pxor xmm1, xmm11            # xmm1 := Y2 xor shuffle(Y2, 78)
    pclmulqdq xmm1, xmm5, 0x00  # xmm1 := H2_Y2_mid

    # - xmm2 := H3_Y3_mid
    pshufd xmm10, xmm6, 78      # xmm10 := shuffle(H3, 78)
    pshufd xmm11, xmm2, 78      # xmm11 := shuffle(Y3, 78)
    pxor xmm6, xmm10            # xmm6 := H3 xor shuffle(H3, 78)
    pxor xmm2, xmm11            # xmm2 := Y3 xor shuffle(Y3, 78)
    pclmulqdq xmm2, xmm6, 0x00  # xmm2 := H3_Y3_mid

    # - xmm3 := H4_Y4_mid
    pshufd xmm10, xmm7, 78      # xmm10 := shuffle(H4, 78)
    pshufd xmm11, xmm3, 78      # xmm11 := shuffle(Y4, 78)
    pxor xmm7, xmm10            # xmm7 := H4 xor shuffle(H4, 78)
    pxor xmm3, xmm11            # xmm3 := Y4 xor shuffle(Y4, 78)
    pclmulqdq xmm3, xmm7, 0x00  # xmm3 := H4_Y4_mid


    # -- Step 1.4: compute C3 + C0 + mid_1 + mid_2 + mid_3 + mid_4
    # copys of C3 and C0
    movdqu xmm10, xmm8          # xmm10 := C0 copy
    movdqu xmm11, xmm9          # xmm11 := C3 copy

    # - xmm10 := C3 + C0 + mid_1 + mid_2 + mid_3 + mid_4
    pxor xmm10, xmm11           # xmm10 := C0 xor C3
    pxor xmm10, xmm0            # xmm10 := H1_Y1_mid xor C0 xor C3
    pxor xmm10, xmm1            # xmm10 := H2_Y2_mid xor H1_Y1_mid xor C0 xor C3
    pxor xmm10, xmm2            # xmm10 := H3_Y3_mid xor H2_Y2_mid xor H1_Y1_mid xor C0 xor C3
    pxor xmm10, xmm3            # xmm10 := H4_Y4_mid xor H3_Y3_mid xor H2_Y2_mid xor H1_Y1_mid xor C0 xor C3


    # -- Step 1.5: final assemble of [C3:C2:C1:C0]
    # copy of xmm10
    movdqu xmm11, xmm10         # xmm11 := xmm10 copy = MID_TERM copy
    pslldq xmm11, 8             # xmm11 := MID_TERM << 64
    psrldq xmm10, 8             # xmm10 := MID_TERM >> 64
    pxor xmm8, xmm11            # xmm8 := C0 xor (MID_TERM << 64)
    pxor xmm9, xmm10            # xmm9 := C3 xor (MID_TERM >> 64)

    # --> xmm8 := [C1:C0] = lo, xmm9 := [C3:C2] = hi
    movdqu xmm1, xmm8           # xmm1 := xmm8 = lo
    movdqu xmm2, xmm1           # xmm2 := xmm1 = xmm8 = lo
    movdqu xmm0, xmm9           # xmm0 := xmm9 = hi
    movdqu xmm5, xmm0           # xmm5 := xmm0 = xmm9 = hi


    # --- Step 2.1: Bitshift <<_epi256 1
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

    # --> Now xmm5 := hi; xmm2 := lo
    movdqa xmm0, xmm2           # xmm0 := xmm2 = copy(lo)
    movdqa xmm1, xmm2           # xmm1 := xmm2 = copy(lo)

    # -- Step 2.2: Load constants
    movdqu xmm6, [rip + Q_r_128_asm]     # xmm6 := Q_r_128
    pxor xmm8, xmm8                      # xmm8 := ZERO (storing ZERO as constant and loading twice is more load than just doing one xor every time)

    # -- Step 2.3: Reduce
    pclmulqdq xmm1, xmm6, 0x00          # xmm1 := x = lo[0] * Q_r[0]
    movdqa xmm7, xmm1                   # xmm7 := xmm1 = copy(x)
    psrldq xmm7, 8                      # xmm7 := x >>_epi128 64
    pxor xmm5, xmm7                     # xmm5 := hi = hi xor (x >>_epi128 64)
    punpcklqdq xmm0, xmm8               # xmm0 := _mm_unpacklo_epi64(lo, ZERO)
    pxor xmm5, xmm0                     # xmm5 := hi = hi xor (_mm512_unpacklo_epi64(lo, ZERO))
    pslldq xmm1, 8                      # xmm1 := x <<_epi128 64
    pxor xmm2, xmm1                     # xmm2 := lo = lo xor (x <<_epi128 64)

    # -- Step 2.4: Reduce
    pclmulqdq xmm6, xmm2, 0x10          # xmm6 := x = lo[1] * Q_r[0]
    pxor xmm5, xmm6                     # xmm5 := hi = hi xor x
    punpckhqdq xmm8, xmm2               # xmm8 := _mm_unpackhi_epi64(ZERO, lo)
    pxor xmm5, xmm8                     # xmm5 := hi = hi xor (_mm512_unpackhi_epi64(ZERO, lo))


    # --- Step 3.0: Store result
    mov rcx, [rsp + 72]                 # rcx := address of __m128i out
    movdqu [rcx], xmm5
    ret


#ifdef __ELF__
.size gfmul_reflected_avx128, .-gfmul_reflected_avx128
.size gfmul_reflected_avx128_times_four, .-gfmul_reflected_avx128_times_four
#endif

.att_syntax
