#include <gcm/gcm.h>
#include <QDebug>
#include <QString>

#include <cpuid.h>
#include <iostream>

bool hasSSSE3() {
    unsigned int eax, ebx, ecx, edx;
    if (__get_cpuid(1, &eax, &ebx, &ecx, &edx)) {
        return (ecx & bit_SSSE3);
    }
    return false;
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
