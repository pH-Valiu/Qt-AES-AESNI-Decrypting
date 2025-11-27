#include <gcm/gcm.h>
#include <aesni/aesni-key-init.h>
#include <qaesencryption.h>
#include <QDebug>
#include <QString>

void singleAESBlock(const unsigned char* in, unsigned char* out, const unsigned char* key, int number_of_rounds){
    __m128i tmp = _mm_loadu_si128((__m128i*) in);
    tmp = _mm_xor_si128(tmp, ((__m128i*)key)[0]);
    for(int i=1; i<number_of_rounds; i++){
        tmp = _mm_aesenc_si128(tmp, ((__m128i*)key)[i]);
    }
    tmp = _mm_aesenclast_si128(tmp, ((__m128i*)key)[number_of_rounds]);
    _mm_storeu_si128((__m128i*)out, tmp);
}

struct GCM_OUT encrypt(const QByteArray &key, const QByteArray &iv, const QByteArray &aad, const QByteArray &p){
    // 0.1: Check instruction set support
    if(!(hasAES() && hasAVX2() && hasPCLMUL())) {
        throw std::runtime_error("CPU does not support required AES-GCM instructions (AES-NI, PCLMUL, AVX2");
    }

    // 0.2: Check max lengths of paramteres following NIST specification (NIST 800-38d)
    if(key.length() != gcm_keyLength || iv.isEmpty() || iv.length() > MAX_IV_LEN || aad.length() > MAX_AAD_LEN || p.length() > MAX_PLAIN_LEN) {
        return GCM_OUT();
    }

    // 1. Key Expansion
    AES_KEY aesKey;
    AES_set_encrypt_key((unsigned char*) key.constData(), gcm_keyLengthBits, &aesKey);
    QByteArray expKey = QByteArray::fromRawData((const char*) aesKey.KEY, gcm_expKeyLength);    // bc expKey is fromRawData, we must not delete aesKey, because expKey contains the data pointer of aesKey. but killing expKey will not kill the key inside

    // 2. IV Expansion



    return GCM_OUT();

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
