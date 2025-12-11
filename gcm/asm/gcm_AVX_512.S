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

.global BSWAP_MASK_asm
.global BSWAP_MASK_TIMES_FOUR_asm
.global BSWAP_MASK_EPI64_asm
.global BSWAP_MASK_EPI64_TIMES_FOUR_asm
.global ONE_asm
.global FOUR_asm
.global ONE_TWO_THREE_FOUR_asm
.global FOUR_TIMES_FOUR_asm

#ifdef __ELF__
.type BSWAP_MASK_asm, @object
.type BSWAP_MASK_EPI64_asm, @object
.type ONE_asm, @object
.type FOUR_asm, @object
#endif

BSWAP_MASK_asm:
    .quad 0x08090a0b0c0d0e0f
    .quad 0x0001020304050607
BSWAP_MASK_TIMES_FOUR_asm:
    .quad 0x08090a0b0c0d0e0f
    .quad 0x0001020304050607
    .quad 0x08090a0b0c0d0e0f
    .quad 0x0001020304050607
    .quad 0x08090a0b0c0d0e0f
    .quad 0x0001020304050607
    .quad 0x08090a0b0c0d0e0f
    .quad 0x0001020304050607
BSWAP_MASK_EPI64_asm:
    .quad 0x0001020304050607
    .quad 0x08090a0b0c0d0e0f
BSWAP_MASK_EPI64_TIMES_FOUR_asm:
    .quad 0x0001020304050607
    .quad 0x08090a0b0c0d0e0f
    .quad 0x0001020304050607
    .quad 0x08090a0b0c0d0e0f
    .quad 0x0001020304050607
    .quad 0x08090a0b0c0d0e0f
    .quad 0x0001020304050607
    .quad 0x08090a0b0c0d0e0f
ONE_asm:
    .quad 0x0000000000000000
    .quad 0x0000000000000001
FOUR_asm:
    .quad 0x0000000000000000
    .quad 0x0000000000000004
ONE_TWO_THREE_FOUR_asm:     # these should correspond to [511:0] =: [0:0x00000004:0:0x00000003:0:0x00000002:0:0x00000001]
    .quad 0x0000000000000000    # these are bits [63-0]
    .quad 0x0000000000000001    # [127-64]
    .quad 0x0000000000000000    # [191-128]
    .quad 0x0000000000000002    # [255-192]
    .quad 0x0000000000000000    # ...
    .quad 0x0000000000000003
    .quad 0x0000000000000000
    .quad 0x0000000000000004
FOUR_TIMES_FOUR_asm:
    .quad 0x0000000000000000
    .quad 0x0000000000000004
    .quad 0x0000000000000000
    .quad 0x0000000000000004
    .quad 0x0000000000000000
    .quad 0x0000000000000004
    .quad 0x0000000000000000
    .quad 0x0000000000000004


#ifdef __ELF__
.size BSWAP_MASK_asm, 16   # tell the linker its size
.size BSWAP_MASK_EPI64_asm, 16
#endif


#--------------------
# ALGORITHM
#--------------------

.text
.intel_syntax noprefix
.global gcm_encrypt_ghash_avx512

#ifdef __ELF__
.type gcm_encrypt_ghash_avx512, @function
#endif

# X_out_ptr         <____ rbp + 104
# H4_ptr            <____ rbp + 96
# H3_ptr            <____ rbp + 88
# H2_ptr            <____ rbp + 80
# H1_ptr            <____ rbp + 72
# X_prev_ptr        <____ rbp + 64
# c_link_ptr        <____ rbp + 56
# p.constData()_ptr <____ rbp + 48
# shadow space r9   <____ rbp + 40
# shadow space r8   <____ rbp + 32
# shadow space rdx  <____ rbp + 24
# shadow space rcx  <____ rbp + 16
# nextInstruction   <____ rbp + 8
# old_rbp           <____ rbp
# stack space [H4, H4, H4, H4]              <____ rbp-64
# stack space [H1, H2, H3, H4]              <____ rbp-128
# stack space [tmp4, tmp3, tmp2, tmp1]      <____ rbp-192
# stack sapce [res4, res3, res2, res1]      <____ rbp-256
# caller saved r12
# caller saved r13
# caller saved r14
# caller saved r15                          <____ rsp
#----------     this is created and destroyed many times inside the method
# caller saved rcx
# caller saved rdx
# caller saved r8
# 512 = 8*16 bytes for caller saved zmm registers for gfmul_reflected_avx512_512
#----------

gcm_encrypt_ghash_avx512:
    # Paramteres needed:
    # - Y0                                  (ptr to value)          (rcx)
    # - key_ptr                             (ptr to array)          (rdx)
    # - aesKey.nr                           (value)                 (r8)
    # - block count (p.length() / 16 / 4)   (value)                 (r9)
    # - p.constData()                       (ptr to array)          (rbp + 48)
    # - c_link                              (ptr to array)          (rbp + 56)
    # - X_prev                              (ptr to value)          (rbp + 64)
    # - H1                                  (ptr to value)          (rbp + 72)
    # - H2                                  (ptr to value)          (rbp + 80)
    # - H3                                  (ptr to value)          (rbp + 88)
    # - H4                                  (ptr to value)          (rbp + 96)
    # - X_out                               (ptr to value)          (rbp + 104)

    # prologue
    push rbp        # this makes rsp-8 and thus again 16 byte aligned, because stack is already 16 byte aligned just before the call instruction which then pushed the rip onto the stack (rsp-8)
    mov rbp, rsp
    sub rsp, 256     # space for 4x 16Byte variables (make sure that rsp is 16 byte aligned - which it is at this point)
    push r12
    push r13
    push r14
    push r15


    # main function

    # static variable loading and ctr1-4 initialization
    vmovdqu xmm0, [rcx]                         # xmm0 := Y0
    vmovdqu xmm1, [rip+ BSWAP_MASK_EPI64_asm]   # xmm1 := BSWAP_MASK_EPI64_asm
    pshufb xmm0, xmm1                           # xmm0 := _mm_shuffle_epi8(Y0, BSWAP_MASK_EPI64_asm) = ctr0
    vinserti32x4 zmm0, zmm0, xmm0, 0            # zmm0 := [RAW, RAW, RAW, xmm0] = [RAW, RAW, RAW, RAW, ctr0]
    vshufi32x4 zmm0, zmm0, zmm0, 0              # zmm0 := [xmm0, xmm0, xmm0, xmm0] = [ctr0, ctr0, ctr0, ctr0]
    vmovdqu64 zmm3, [rip + ONE_TWO_THREE_FOUR_asm]# zmm3 := [FOUR, THREE, TWO, ONE]
    vpaddd zmm0, zmm0, zmm3                     # zmm0 := [xmm0 + 4, xmm0 + 3, xmm0 + 2, xmm0 + 1] = [ctr4, ctr3, ctr2, ctr1]

    vmovdqu64 zmm3, [rip + BSWAP_MASK_EPI64_TIMES_FOUR_asm]   # zmm3 := BSWAP_MASK_EPI64_TIMES_FOUR_asm
    vmovdqu64 zmm4, [rip + FOUR_TIMES_FOUR_asm]   # zmm4 := [FOUR, FOUR, FOUR, FOUR]
    xor r10, r10                                # r10 := i = 0      (block loop variable)
    shl r9, 6                                   # r9 := block count * 16 * 4 = block count * 64 =: byte count (for all the full blocks)

    cmp r10, r9
    je .L_epilogue                              # early exit in case that we do not have a single 4times block of data

    # otherwise: preload the aes keys (zmm15 - zmm23) = KEY[1] - KEY[9] and zmm5 = KEY[0]
    mov r13, rdx                                # r13 := rdx = &key_ptr[0] =: j
    vmovdqu xmm5, [r13]                         # xmm5 := key_ptr[0]
    vinserti32x4 zmm5, zmm5, xmm5, 0            # zmm5 := [RAW, RAW, RAW, xmm5]
    vshufi32x4 zmm5, zmm5, zmm5, 0              # zmm5 := [xmm5, xmm5, xmm5, xmm5] = [key_ptr[0], key_ptr[0], key_ptr[0], key_ptr[0]]

    add r13, 16
    vmovdqu64 xmm15, [r13]                        # xmm15 := key_ptr[16] (key_ptr is byte pointer) = KEY[1]
    vinserti32x4 zmm15, zmm15, xmm15, 0         # zmm15 := [RAW, RAW, RAW, xmm15]
    vshufi32x4 zmm15, zmm15, zmm15, 0           # zmm15 := [xmm15, xmm15, xmm15, xmm15] = [key_ptr[16], key_ptr[16], key_ptr[16], key_ptr[16]]

    add r13, 16
    vmovdqu64 xmm16, [r13]                        # xmm16 := key_ptr[32] = KEY[2]
    vinserti32x4 zmm16, zmm16, xmm16, 0
    vshufi32x4 zmm16, zmm16, zmm16, 0

    add r13, 16
    vmovdqu64 xmm17, [r13]                        # xmm17 := key_ptr[48] = KEY[3]
    vinserti32x4 zmm17, zmm17, xmm17, 0
    vshufi32x4 zmm17, zmm17, zmm17, 0

    add r13, 16
    vmovdqu64 xmm18, [r13]                        # xmm18 := key_ptr[64] = KEY[4]
    vinserti32x4 zmm18, zmm18, xmm18, 0
    vshufi32x4 zmm18, zmm18, zmm18, 0

    add r13, 16
    vmovdqu64 xmm19, [r13]                        # xmm19 := key_ptr[80] = KEY[5]
    vinserti32x4 zmm19, zmm19, xmm19, 0
    vshufi32x4 zmm19, zmm19, zmm19, 0

    add r13, 16
    vmovdqu64 xmm20, [r13]                        # xmm20 := key_ptr[96] = KEY[6]
    vinserti32x4 zmm20, zmm20, xmm20, 0
    vshufi32x4 zmm20, zmm20, zmm20, 0

    add r13, 16
    vmovdqu64 xmm21, [r13]                        # xmm21 := key_ptr[112] = KEY[7]
    vinserti32x4 zmm21, zmm21, xmm21, 0
    vshufi32x4 zmm21, zmm21, zmm21, 0

    add r13, 16
    vmovdqu64 xmm22, [r13]                        # xmm22 := key_ptr[128] = KEY[8]
    vinserti32x4 zmm22, zmm22, xmm22, 0
    vshufi32x4 zmm22, zmm22, zmm22, 0

    add r13, 16
    vmovdqu64 xmm23, [r13]                        # xmm23 := key_ptr[144] = KEY[9]
    vinserti32x4 zmm23, zmm23, xmm23, 0
    vshufi32x4 zmm23, zmm23, zmm23, 0

    # other necessary variable for last aes loops
    sub r8, 1                                   # r8 := aesKey.nr - 1
    shl r8, 4                                   # r8 := (aesKey.nr - 1) * 16

    vmovdqu64 zmm10, [rip + BSWAP_MASK_TIMES_FOUR_asm]  # zmm10 := BSWAP_MASK_TIMES_FOUR_asm
    vpxord zmm11, zmm11, zmm11
    vmovdqu xmm6, [rbp + 72]                    # xmm6 := H1
    vmovdqu xmm7, [rbp + 80]                    # xmm7 := H2
    vmovdqu xmm8, [rbp + 88]                    # xmm8 := H3
    vmovdqu xmm9, [rbp + 96]                    # xmm9 := H4
    vinserti32x4 zmm25, zmm25, xmm9, 0          # zmm25 := [RAW, RAW, RAW, H4]
    vshufi32x4 zmm25, zmm25, zmm24, 0           # zmm25 := [H4, H4, H4, H4]
    vmovdqu64 [rbp - 64], zmm25                 # store [H4, H4, H4, H4] into [rbp - 64]
    vinserti32x4 zmm26, zmm26, xmm9, 0
    vinserti32x4 zmm26, zmm26, xmm8, 1
    vinserti32x4 zmm26, zmm26, xmm7, 2
    vinserti32x4 zmm26, zmm26, xmm6, 3          # zmm26 := [H1, H2, H3, H4]
    vmovdqu64 [rbp - 128], zmm26               # store [H1, H2, H3, H4] into [rbp - 128]

    # start main block loop
    jmp .L_block_loop

    # we state that zmm11 shall hold the res4-1 registers
    # zmm0 holds ctr4-1
    # zmm1 holds tmp4-1
    # zmm3 holds BSWAP_MASK_EPI64_TIMES_FOUR_asm
    # zmm4 holds FOUR_TIMES_FOUR_asm
    # zmm5 holds [KEY[0], KEY[0], KEY[0], KEY[0]]
    # zmm6 free
    # zmm7 free
    # zmm8 free
    # zmm9 free
    # zmm10 holds BSWAP_MASK_TIMES_FOUR_asm
    # zmm11 holds res4-1
    # zmm12 free
    # zmm13 free
    # zmm14 free
    # zmm15 holds KEY[1]
    # zmm16 holds KEY[2]
    # zmm17 holds KEY[3]
    # zmm18 holds KEY[4]
    # zmm19 holds KEY[5]
    # zmm20 holds KEY[6]
    # zmm21 holds KEY[7]
    # zmm22 holds KEY[8]
    # zmm23 holds KEY[9]
    # zmm24 shall hold KEY[10-aesKey.nr]
    # zmm25 holds [H4, H4, H4, H4]  (also present at [rbp - 64])
    # zmm26 holds [H1, H2, H3, H4]  (also present at [rbp - 128])
    #
    # r8 holds (aesKey.nr - 1) * 16
    # r9 holds (p.length() / 16 / 4) * 16 * 4
    # r10 holds i
    # r11 free
    # r12 free
    # r13 holds j (free)
    # r14
    # r15

.L_block_loop:
    cmp r10, r9                                 # if r10 := i == block_count * 16 = p.lenght()/16/4 * 16 =: r9
    je .L_epilogue                              # jump to epilogue if equal
    vpshufb zmm1, zmm0, zmm3                    # zmm1 := [shuffle(xmm0+4), shuffle(xmm0+3), shuffle(xmm0+2), shuffle(xmm0+1)] = [shuffle(ctr4), shuffle(ctr3), shuffle(ctr2), shuffle(ctr1)] = [tmp4, tmp3, tmp2, tmp1]
    vpaddd zmm0, zmm0, zmm4                     # zmm0 := [xmm0+4 + 4, xmm0+3 + 4, xmm0+2 + 4, xmm0+1 + 4] = [ctr4 + 4, ctr3 + 4, ctr2 + 4, ctr1 + 4]

    # reminder: zmm1 := [tmp4, tmp3, tmp2, tmp1]
    # apply round key KEY[0] =: zmm5
    vpxord zmm1, zmm1, zmm5                     # zmm1 := zmm1 xor zmm5 = [tmp4, tmp3, tmp2, tmp1] xor [key_ptr[0], key_ptr[0], key_ptr[0], key_ptr[0]]

    # apply awes rounds 1 - 9
    vaesenc zmm1, zmm1, zmm15                         # zmm1 := AES(zmm1, zmm15) = AES(zmm1, KEY[1])
    vaesenc zmm1, zmm1, zmm16
    vaesenc zmm1, zmm1, zmm17
    vaesenc zmm1, zmm1, zmm18
    vaesenc zmm1, zmm1, zmm19
    vaesenc zmm1, zmm1, zmm20
    vaesenc zmm1, zmm1, zmm21
    vaesenc zmm1, zmm1, zmm22
    vaesenc zmm1, zmm1, zmm23                         # zmm1 := AES(zmm1, zmm23) = AES(zmm1, KEY[9])

    # check here if more rounds have to be performed (e.g. AES-196 and AES-256)
    mov r13, 144                                # r13 := 144 = 9 * 16 (144 is the last aes key before the last last key for AES-128)
    jmp .L_aes_additional_loops

.L_aes_additional_loops:
    cmp r13, r8
    je .L_aes_last_round

    # additional aes loops
    add r13, 16
    vmovdqu64 xmm24, [r13 + rdx]                    # xmm24 := key_ptr[144 + 16*(last_rounds_index+1)] = KEY[9 + (last_rounds_index+1)]
    vinserti32x4 zmm24, zmm24, xmm24, 0
    vshufi32x4 zmm24, zmm24, zmm24, 0
    vaesenc zmm1, zmm1, zmm24                         # zmm1 := AES(zmm1, zmm24) = AES(zmm1, KEY[10 + ...])

    jmp .L_aes_additional_loops


.L_aes_last_round:
    add r13, 16
    vmovdqu64 xmm24, [r13 + rdx]                  # xmm24 := key_ptr[160] = KEY[10]
    vinserti32x4 zmm24, zmm24, xmm24, 0
    vshufi32x4 zmm24, zmm24, zmm24, 0
    vaesenclast zmm1, zmm1, zmm24                     # zmm1 := AES_LAST(zmm1, zmm24) = AES_LAST(zmm1, KEY[last])  (should be 10 for AES-128, 12 for AES-196, 14 for AES-256)


    # load data
    mov r11, [rbp + 48]                         # r11 := p.constData()
    #lea r12, [r10 + 0]                          # r12 := r10 + 0 = i
    #vmovdqu xmm6, [r11 + r12]                   # xmm6 := ((__m128i*)p.constData())[i]
    #lea r12, [r10 + 16]                         # r12 := r10 + 16 = i+1
    #vmovdqu xmm7, [r11 + r12]                   # xmm7 := ((__m128i*)p.constData())[i+1]
    #lea r12, [r10 + 32]                         # r12 := r10 + 32 = i+2
    #vmovdqu xmm8, [r11 + r12]                   # xmm8 := ((__m128i*)p.constData())[i+2]
    #lea r12, [r10 + 48]                         # r12 := r10 + 48 = i+3
    #vmovdqu xmm9, [r11 + r12]                   # xmm9 := ((__m128i*)p.constData())[i+3]

    # merge data to one zmm register (zmm6)
    #vinserti64x2 zmm6, zmm6, xmm6, 0
    #vinserti64x2 zmm6, zmm6, xmm7, 1
    #vinserti64x2 zmm6, zmm6, xmm8, 2
    #vinserti64x2 zmm6, zmm6, zmm9, 3            # zmm6 := [xmm9, xmm8, xmm7, xmm6]

    vmovdqu64 zmm6, [r11 + r10]                  # zmm6 = [p[i+3], p[i+2], p[i+1], p[i]]

    # encrypt data with zmm1 being the AES result
    vpxord zmm1, zmm1, zmm6                     # zmm1 := zmm1 xor zmm6 = [tmp4, tmp3, tmp2, tmp1] xor [p[i+3], p[i+2], p[i+1], p[i]] (this is the actual encryption of the data)

    # store encrypted data in c_link
    mov r14, [rbp + 56]
    vmovdqu64 [r14 + r10], zmm1                 # store zmm1 being the encrypted cipher into stack storage at c_link ptr

    # shuffle tmp4-1 for GHASH
    vpshufb zmm6, zmm6, zmm10                   # zmm6 := shuffle(zmm6, BSWAP_MASK) = [shuffle_epi128(tmp4), shuffle_epi128(tmp3), shuffle_epi128(tmp2), shuffle_epi128(tmp1)]

    xor r15, r15
    cmp r10, r15
    je .L_single_xor_X_prev
    vpxord zmm6, zmm6, zmm11                    # zmm6 := zmm6 xor zmm11 = tmp4-1 xor res4-1
    jmp .L_after_xor_X_prev_or_res

.L_single_xor_X_prev:
    vpxord zmm12, zmm12, zmm12
    vmovdqu xmm13, [rbp + 64]
    vinserti64x2 zmm12, zmm12, xmm13, 0
    vpxord zmm6, zmm6, zmm12

    jmp .L_after_xor_X_prev_or_res

.L_after_xor_X_prev_or_res:
    mov r11, r10                                # r11 := r10 = i
    add r11, 64                                 # r11 := r11 + 64
    cmp r11, r9                                 # r11 ?= r9 : i+1 ?= p_byte_length
    je .L_last_block_round

    # for the call, we need rcx, rdx and r8 to be filled
    # rcx = a = ptr to [tmp4, tmp3, tmp2, tmp1]
    # rdx = b = ptr to [H4, H4, H4, H4]
    # r8 = res = ptr to [res4, res3, res2, res1]

    # store rcx, rdx, r8 and restore after method call
    push rcx
    push rdx
    push r8
    push r9         # just so that we have 16 byte alignment again

    # the other r.. registers which are caller saved do not need to be saved as gfmul_reflected_avx512_512 does not use them

    vmovdqu64 [rbp - 192], zmm1      # store [tmp4, tmp3, tmp2, tmp1] into [rbp-192]
    lea rcx, [rbp - 192]            # rcx := ptr to tmp4-1
    lea rdx, [rbp - 64]             # rdx := value of rbp - 64 = ptr to [H4, H4, H4, H4]
    vmovdqu64 [rbp - 256], zmm11     # store [res4, res3, res2, res1] into [rbp-256]
    lea r8, [rbp - 256]             # t8 := ptr to res4-1

    # store zmm registers that gfmul_reflected_avx512_512 uses
    sub rsp, 512
    vmovdqu64 [rsp], zmm0
    vmovdqu64 [rsp+64], zmm1
    vmovdqu64 [rsp+128], zmm2
    vmovdqu64 [rsp+192], zmm3
    vmovdqu64 [rsp+256], zmm4
    vmovdqu64 [rsp+320], zmm5
    vmovdqu64 [rsp+384], zmm6
    vmovdqu64 [rsp+448], zmm7
    call gfmul_reflected_avx512_512
    vmovdqu64 zmm7, [rsp+448]
    vmovdqu64 zmm6, [rsp+384]
    vmovdqu64 zmm5, [rsp+320]
    vmovdqu64 zmm4, [rsp+256]
    vmovdqu64 zmm3, [rsp+192]
    vmovdqu64 zmm2, [rsp+128]
    vmovdqu64 zmm1, [rsp+64]
    vmovdqu64 zmm0, [rsp]
    add rsp, 512

    # load back registers rcx, rdx, r8, r9
    pop r9
    pop r8
    pop rdx
    pop rcx

    # store result back into zmm11
    vmovdqu64 zmm11, [rbp - 256]                # zmm11 := load from stack = [res4, res3, res2, res1]


    # increment loop index
    add r10, 64                                 # i++ (+64 = +16*4)
    jmp .L_block_loop                           # loop back


.L_last_block_round:
    # for the call, we need rcx, rdx and r8 to be filled
    # rcx = a = ptr to [tmp4, tmp3, tmp2, tmp1]
    # rdx = b = ptr to [H1, H2, H3, H4]
    # r8 = res = ptr to [res4, res3, res2, res1]

    # we need not to store register r..., as gfmul_reflected_avx512_512 does not use them and this is leave call.
    vmovdqu64 [rbp - 192], zmm1      # store [tmp4, tmp3, tmp2, tmp1] into [rbp-192]
    lea rcx, [rbp - 192]            # rcx := ptr to tmp4-1
    lea rdx, [rbp - 128]            # rdx := value of rbp - 128 = ptr to [H1, H2, H3, H4]
    vmovdqu64 [rbp - 256], zmm11     # store [res4, res3, res2, res1] into [rbp-256]
    lea r8, [rbp - 256]             # t8 := ptr to res4-1

    # store zmm registers that gfmul_reflected_avx512_512 uses
    sub rsp, 512
    vmovdqu64 [rsp], zmm0
    vmovdqu64 [rsp+64], zmm1
    vmovdqu64 [rsp+128], zmm2
    vmovdqu64 [rsp+192], zmm3
    vmovdqu64 [rsp+256], zmm4
    vmovdqu64 [rsp+320], zmm5
    vmovdqu64 [rsp+384], zmm6
    vmovdqu64 [rsp+448], zmm7
    call gfmul_reflected_avx512_512
    vmovdqu64 zmm7, [rsp+448]
    vmovdqu64 zmm6, [rsp+384]
    vmovdqu64 zmm5, [rsp+320]
    vmovdqu64 zmm4, [rsp+256]
    vmovdqu64 zmm3, [rsp+192]
    vmovdqu64 zmm2, [rsp+128]
    vmovdqu64 zmm1, [rsp+64]
    vmovdqu64 zmm0, [rsp]
    add rsp, 512

    vmovdqu64 zmm11, [r8]                       # zmm11 := res4-1 (as r8 holds address of stack storage for res4-1 given to gfmul_reflected method)
    vextracti64x2 xmm6, zmm11, 0                # xmm6 := res1
    vextracti64x2 xmm7, zmm11, 1                # xmm7 := res2
    vextracti64x2 xmm8, zmm11, 2                # xmm8 := res3
    vextracti64x2 xmm9, zmm11, 3                # xmm9 := res4
    pxor xmm6, xmm7
    pxor xmm6, xmm8
    pxor xmm6, xmm9                             # xmm6 := xmm6 xor xmm7 xor xmm8 xor xmm9
    mov r8, [rbp + 104]                         # r8 := X_out_ptr
    vmovdqu [r8], xmm6                          # store xmm6 := X in [r8] = X_out
    jmp .L_epilogue


.L_epilogue:
    # epilogue
    pop r15
    pop r14
    pop r13
    pop r12
    mov rsp, rbp
    pop rbp
    ret


#ifdef __ELF__
.size gcm_encrypt_ghash_avx512, .-gcm_encrypt_ghash_avx512
#endif

.att_syntax
