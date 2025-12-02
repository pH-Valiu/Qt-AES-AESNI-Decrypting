#include <gcm/gcm.h>
#include <gcm/gfmul.h>
#include <gcm/gcm_compatability.h>
#include <aesni/aesni-key-init.h>
#include <gcm/benchmarkutil.h>

#include <smmintrin.h>
#include <qaesencryption.h>
#include <QDebug>
#include <QString>
static inline void assert_m128(__m128i a, __m128i a_assert, const QString& name)
{
    __m128i diff = _mm_xor_si128(a, a_assert);
    __m128i mask = _mm_set1_epi32(0xFFFFFFFF);

    if (_mm_test_all_zeros(mask, diff)) {
        qInfo() << "Assertion (" << name << "): OK";
    } else {
        qWarning() << "Assertion (" << name << "): FAILED";
        qWarning() << "Expected: "<<print128_hex_lanes(a_assert);
        qWarning() << "Actual: "<<print128_hex_lanes(a);
    }
}

void singleAESBlock(const unsigned char* in, unsigned char* out, const unsigned char* key, int number_of_rounds){
    __m128i tmp = _mm_loadu_si128((__m128i*) in);   // take first 16 bytes from in and store in tmp pointer
    tmp = _mm_xor_si128(tmp, ((__m128i*)key)[0]);   // apply Round 0 key
    for(int i=1; i<number_of_rounds-1; i+=2){
        tmp = _mm_aesenc_si128(tmp, ((__m128i*)key)[i]);    // Rounds: [1 - (n-1)]
        tmp = _mm_aesenc_si128(tmp, ((__m128i*)key)[i+1]);    // Rounds: [1 - (n-1)]
    }
    tmp = _mm_aesenc_si128(tmp, ((__m128i*)key)[number_of_rounds-1]);    // Rounds: [1 - (n-1)]
    tmp = _mm_aesenclast_si128(tmp, ((__m128i*)key)[number_of_rounds]); // Round: n
    _mm_storeu_si128((__m128i*)out, tmp);           // store tmp in first 16 bytes of out pointer
}

__m128i singleAESBlock(const __m128i& in, const __m128i* const key, int number_of_rounds){
    __m128i out = _mm_xor_si128(in, key[0]);
    for(int i=1; i<number_of_rounds-1; i+=2){
        out = _mm_aesenc_si128(out, key[i]);
        out = _mm_aesenc_si128(out, key[i+1]);
    }
    out = _mm_aesenc_si128(out, key[number_of_rounds-1]);
    return _mm_aesenclast_si128(out, key[number_of_rounds]);
}

__m128i singleAESBlock(const __m128i& in, const AES_KEY& key){
    __m128i* keys = (__m128i*) key.KEY;
    __m128i out = _mm_xor_si128(in, keys[0]);
    for(int i=1; i<key.nr-1; i+=2){
        out = _mm_aesenc_si128(out, keys[i]);
        out = _mm_aesenc_si128(out, keys[i+1]);
    }
    out = _mm_aesenc_si128(out, keys[key.nr-1]);
    return _mm_aesenclast_si128(out, keys[key.nr]);
}

GCM_OUT encrypt(const QByteArray &key, const QByteArray &iv, const QByteArray &aad, const QByteArray &p){
    // 0.1: Check instruction set support
    if(!(AES_GCM_Compatability::hasAES() && AES_GCM_Compatability::hasAVX2() && AES_GCM_Compatability::hasPCLMUL())) {
        throw std::runtime_error("CPU does not support required AES-GCM instructions (AES-NI, PCLMUL, AVX2");
    }

    // 0.2: Check max lengths of paramteres following NIST specification (NIST 800-38d)
    if(key.length() != gcm_keyLength || iv.isEmpty() || iv.length() > MAX_IV_LEN || aad.length() > MAX_AAD_LEN || p.length() > MAX_PLAIN_LEN) {
        return GCM_OUT();
    }

    // 0.3: Key Expansion
    AES_KEY aesKey;
    AES_set_encrypt_key((unsigned char*) key.constData(), gcm_keyLengthBits, &aesKey);

    // 0.4 Variable declarations
    __m128i Y0, tmp1, H, T, ctr1;
    __m128i last_block = ZERO;
    __m128i X = ZERO;           // X is used for the running GHASH calculation state
    unsigned int i, j;

    GCM_OUT out;
    out.c.resize(p.length());
    out.t.resize(16);
    char* c_link = out.c.data();

    // 1: H = AES_k(0^128)
    H = singleAESBlock(ZERO, aesKey);
    H = _mm_shuffle_epi8(H, BSWAP_MASK);

    // 2. IV Expansion
    if(iv.length() == 12){
        Y0 = _mm_loadu_si128((__m128i*) iv.constData());
        Y0 = _mm_insert_epi32(Y0, 0x01000000, 3);       // this is such that in memory a0:| iv0, iv1, iv2, iv3, iv4, ..., iv11, iv12, 00, 00, 00, 01 |:a15
    } else {        // We have to apply GHASH(IV||0^{remainingBytesForFullBlock}||0^32||iv.bitlength())
        Y0 = ZERO;
        for(i=0; i < iv.length()/16; i++){  // We do first all full 16 byte blocks
            tmp1 = _mm_loadu_si128(&((__m128i*)iv.constData())[i]); // the _mm_loadu_si128 is necessary because iv might not be 16Byte aligned
            tmp1 = _mm_shuffle_epi8(tmp1, BSWAP_MASK);
            Y0 = _mm_xor_si128(Y0, tmp1);
            Y0 = gfmul_reflected(Y0, H);
        }
        if(iv.length() % 16){
            for(j=0; j < iv.length() % 16; j++)
                ((unsigned char*)&last_block)[j] = iv[i*16+j];
            tmp1 = last_block;
            tmp1 = _mm_shuffle_epi8(tmp1, BSWAP_MASK);
            Y0 = _mm_xor_si128(Y0, tmp1);
            Y0 = gfmul_reflected(Y0, H);
        }
        tmp1 = _mm_insert_epi64(tmp1, iv.length()*8, 0);
        tmp1 = _mm_insert_epi64(tmp1, 0, 1);
        Y0 = _mm_xor_si128(Y0, tmp1);
        Y0 = gfmul_reflected(Y0, H);
        Y0 = _mm_shuffle_epi8(Y0, BSWAP_MASK);      // this brings it back into "normal" world
    }

    // 3. T_pre computation
    T = singleAESBlock(Y0, aesKey);     // this will be xor-ed onto the full GHASH output to form the final tag

    // 4. GHASH(aad)
    for(i=0; i<aad.length()/16; i++){       // first apply GHASH on all full blocks
        tmp1 = _mm_loadu_si128(&((__m128i*)aad.constData())[i]);    // the _mm_loadu_si128 is necessary because the data might not be 16Byte aligned
        tmp1 = _mm_shuffle_epi8(tmp1, BSWAP_MASK);
        X = _mm_xor_si128(X, tmp1);
        X = gfmul_reflected(X, H);
    }
    if(aad.length() % 16){                  // apply GHASH on the remaining block if necessary
        last_block = ZERO;
        for(j=0; j<aad.length() % 16; j++){
            ((unsigned char*) &last_block)[j] = aad[i*16 + j];
        }
        tmp1 = _mm_shuffle_epi8(last_block, BSWAP_MASK);
        X = _mm_xor_si128(X, tmp1);
        X = gfmul_reflected(X, H);
    }

    // 5. Ciphertext computation
    ctr1 = Y0;
    for(i=0; i<p.length() / 16; i++) {      // we traverse all **full** 16byte blocks of p
        ctr1 = _mm_shuffle_epi8(ctr1, BSWAP_EPI64_MASK);
        ctr1 = _mm_add_epi32(ctr1, ONE);                    // increase ctr1 <- ctr1 + 1
        ctr1 = _mm_shuffle_epi8(ctr1, BSWAP_EPI64_MASK);

        tmp1 = singleAESBlock(ctr1, aesKey);                // encrypt (j0 + 1) -> tmp1

        // the _mm_loadu_si128 is necessary because the data might not be 16Byte aligned
        tmp1 = _mm_xor_si128(tmp1, _mm_loadu_si128(&((__m128i*)p.constData())[i]));     // xor tmp1 with p block -> c block

        // the _mm_store_si128 is necessary because the destination might not be 16Byte aligned
        _mm_storeu_si128(&((__m128i*)c_link)[i], tmp1);           // store c block in out

        tmp1 = _mm_shuffle_epi8(tmp1, BSWAP_MASK);          // bring c block to reverse world

        X = _mm_xor_si128(X, tmp1);
        X = gfmul_reflected(X, H);                          // update GHASH state with currently computed ciphertext C
    }
    if(p.length() % 16){            // handle last block if necessary
        ctr1 = _mm_shuffle_epi8(ctr1, BSWAP_EPI64_MASK);
        ctr1 = _mm_add_epi32(ctr1, ONE);
        ctr1 = _mm_shuffle_epi8(ctr1, BSWAP_EPI64_MASK);

        tmp1 = singleAESBlock(ctr1, aesKey);     // encrypt (j0 + 1) -> tmp1

        tmp1 = _mm_xor_si128(tmp1, _mm_loadu_si128(&((__m128i*)p.constData())[i]));
        last_block = tmp1;
        for(j=0; j<p.length() % 16; j++){
            c_link[i*16 + j] = ((unsigned char*) &last_block)[j];
        }
        for(j; j<16; j++){
            ((unsigned char*) &last_block)[j] = 0;
        }
        tmp1 = _mm_shuffle_epi8(last_block, BSWAP_MASK);
        X = _mm_xor_si128(X, tmp1);
        X = gfmul_reflected(X, H);
    }

    // 6. Final Tag
    tmp1 = _mm_insert_epi64(tmp1, p.length()*8, 0);
    tmp1 = _mm_insert_epi64(tmp1, aad.length()*8, 1);

    X = _mm_xor_si128(X, tmp1);
    X = gfmul_reflected(X, H);
    X = _mm_shuffle_epi8(X, BSWAP_MASK);        // bring final GHASH state back to normal world
    T = _mm_xor_si128(X, T);
    _mm_storeu_si128((__m128i*)out.t.data(), T);

    return out;
}




GCM_OUT encrypt_times_four(const QByteArray &key, const QByteArray &iv, const QByteArray &aad, const QByteArray &p){
    // 0.1: Check instruction set support
    if(!(AES_GCM_Compatability::hasAES() && AES_GCM_Compatability::hasAVX2() && AES_GCM_Compatability::hasPCLMUL())) {
        throw std::runtime_error("CPU does not support required AES-GCM instructions (AES-NI, PCLMUL, AVX2");
    }

    // 0.2: Check max lengths of paramteres following NIST specification (NIST 800-38d)
    if(key.length() != gcm_keyLength || iv.isEmpty() || iv.length() > MAX_IV_LEN || aad.length() > MAX_AAD_LEN || p.length() > MAX_PLAIN_LEN) {
        return GCM_OUT();
    }

    // 0.3: Key Expansion
    AES_KEY aesKey;
    AES_set_encrypt_key((unsigned char*) key.constData(), gcm_keyLengthBits, &aesKey);
    __m128i* key_ptr = (__m128i*) aesKey.KEY;

    // 0.4 Variable declarations
    __m128i Y0, tmp1, tmp2, tmp3, tmp4, H, T, ctr1, ctr2, ctr3, ctr4;
    __m128i last_block = ZERO;
    __m128i X = ZERO;           // X is used for the running GHASH calculation state
    unsigned int i, j, k;

    GCM_OUT out;
    out.c.resize(p.length());
    out.t.resize(16);
    char* c_link = out.c.data();

    // 1: H = AES_k(0^128) and IV Expansion and T_pre computation
    if(iv.length() == 12){
        Y0 = _mm_loadu_si128((__m128i*) iv.constData());
        Y0 = _mm_insert_epi32(Y0, 0x01000000, 3);       // this is such that in memory a0:| iv0, iv1, iv2, iv3, iv4, ..., iv11, iv12, 00, 00, 00, 01 |:a15

        // encrypt H and Y0
        tmp1 = Y0;
        tmp2 = ZERO;        // future H
        tmp1 = _mm_xor_si128(tmp1, key_ptr[0]);
        tmp2 = _mm_xor_si128(tmp2, key_ptr[0]);
        for(i=1; i < aesKey.nr-1; i+=2){
            tmp1 = _mm_aesenc_si128(tmp1, key_ptr[i]);
            tmp2 = _mm_aesenc_si128(tmp2, key_ptr[i]);
            tmp1 = _mm_aesenc_si128(tmp1, key_ptr[i+1]);
            tmp2 = _mm_aesenc_si128(tmp2, key_ptr[i+1]);
        }
        tmp1 = _mm_aesenc_si128(tmp1, key_ptr[aesKey.nr-1]);
        tmp2 = _mm_aesenc_si128(tmp2, key_ptr[aesKey.nr-1]);
        tmp1 = _mm_aesenclast_si128(tmp1, key_ptr[aesKey.nr]);
        tmp2 = _mm_aesenclast_si128(tmp2, key_ptr[aesKey.nr]);

        H = _mm_shuffle_epi8(tmp2, BSWAP_MASK);
        T = tmp1;

    } else {        // We have to apply GHASH(IV||0^{remainingBytesForFullBlock}||0^32||iv.bitlength())
        // compute H
        tmp2 = _mm_xor_si128(ZERO, key_ptr[0]);
        for(i=1; i<aesKey.nr-1; i+=2){
            tmp2 = _mm_aesenc_si128(tmp2, key_ptr[i]);
            tmp2 = _mm_aesenc_si128(tmp2, key_ptr[i+1]);
        }
        tmp2 = _mm_aesenc_si128(tmp2, key_ptr[aesKey.nr-1]);
        tmp2 = _mm_aesenclast_si128(tmp2, key_ptr[aesKey.nr]);
        H = _mm_shuffle_epi8(tmp2, BSWAP_MASK);

        Y0 = ZERO;
        for(i=0; i < iv.length()/16; i++){  // We do first all full 16 byte blocks
            tmp1 = _mm_loadu_si128(&((__m128i*)iv.constData())[i]);
            tmp1 = _mm_shuffle_epi8(tmp1, BSWAP_MASK);
            Y0 = _mm_xor_si128(Y0, tmp1);
            Y0 = gfmul_reflected(Y0, H);
        }
        if(iv.length() % 16){
            for(j=0; j < iv.length() % 16; j++)
                ((unsigned char*)&last_block)[j] = iv[i*16+j];
            tmp1 = last_block;
            tmp1 = _mm_shuffle_epi8(tmp1, BSWAP_MASK);
            Y0 = _mm_xor_si128(Y0, tmp1);
            Y0 = gfmul_reflected(Y0, H);
        }
        tmp1 = _mm_insert_epi64(tmp1, iv.length()*8, 0);
        tmp1 = _mm_insert_epi64(tmp1, 0, 1);
        Y0 = _mm_xor_si128(Y0, tmp1);
        Y0 = gfmul_reflected(Y0, H);
        Y0 = _mm_shuffle_epi8(Y0, BSWAP_MASK);      // this brings it back into "normal" world

        // compute T_pre
        tmp1 = _mm_xor_si128(Y0, key_ptr[0]);
        for(i=1; i<aesKey.nr-1; i+=2){
            tmp1 = _mm_aesenc_si128(tmp1, key_ptr[i]);
            tmp1 = _mm_aesenc_si128(tmp1, key_ptr[i+1]);
        }
        tmp1 = _mm_aesenc_si128(tmp1, key_ptr[aesKey.nr-1]);
        T = _mm_aesenclast_si128(tmp1, key_ptr[aesKey.nr]);
    }

    // 4. GHASH(aad)
    for(i=0; i<aad.length()/16; i++){       // first apply GHASH on all full blocks
        tmp1 = _mm_loadu_si128(&((__m128i*)aad.constData())[i]);
        tmp1 = _mm_shuffle_epi8(tmp1, BSWAP_MASK);
        X = _mm_xor_si128(X, tmp1);
        X = gfmul_reflected(X, H);
    }
    if(aad.length() % 16){                  // apply GHASH on the remaining block if necessary
        last_block = ZERO;
        for(j=0; j<aad.length() % 16; j++){
            ((unsigned char*) &last_block)[j] = aad[i*16 + j];
        }
        tmp1 = _mm_shuffle_epi8(last_block, BSWAP_MASK);
        X = _mm_xor_si128(X, tmp1);
        X = gfmul_reflected(X, H);
    }

    // 5. Ciphertext computation
    ctr1 = _mm_shuffle_epi8(Y0, BSWAP_EPI64_MASK);
    ctr1 = _mm_add_epi32(ctr1, ONE);
    ctr2 = _mm_add_epi32(ctr1, ONE);
    ctr3 = _mm_add_epi32(ctr2, ONE);
    ctr4 = _mm_add_epi32(ctr3, ONE);
    for(i=0; i<p.length() / 16 / 4; i++) {
        tmp1 = _mm_shuffle_epi8(ctr1, BSWAP_EPI64_MASK);
        tmp2 = _mm_shuffle_epi8(ctr2, BSWAP_EPI64_MASK);
        tmp3 = _mm_shuffle_epi8(ctr3, BSWAP_EPI64_MASK);
        tmp4 = _mm_shuffle_epi8(ctr4, BSWAP_EPI64_MASK);

        ctr1 = _mm_add_epi32(ctr1, FOUR);
        ctr2 = _mm_add_epi32(ctr2, FOUR);
        ctr3 = _mm_add_epi32(ctr3, FOUR);
        ctr4 = _mm_add_epi32(ctr4, FOUR);
        tmp1 =_mm_xor_si128(tmp1, key_ptr[0]);
        tmp2 =_mm_xor_si128(tmp2, key_ptr[0]);
        tmp3 =_mm_xor_si128(tmp3, key_ptr[0]);
        tmp4 =_mm_xor_si128(tmp4, key_ptr[0]);
        for(j=1; j < aesKey.nr-1; j+=2){
            tmp1 = _mm_aesenc_si128(tmp1, key_ptr[j]);
            tmp2 = _mm_aesenc_si128(tmp2, key_ptr[j]);
            tmp3 = _mm_aesenc_si128(tmp3, key_ptr[j]);
            tmp4 = _mm_aesenc_si128(tmp4, key_ptr[j]);
            tmp1 = _mm_aesenc_si128(tmp1, key_ptr[j+1]);
            tmp2 = _mm_aesenc_si128(tmp2, key_ptr[j+1]);
            tmp3 = _mm_aesenc_si128(tmp3, key_ptr[j+1]);
            tmp4 = _mm_aesenc_si128(tmp4, key_ptr[j+1]);
        }
        tmp1 = _mm_aesenc_si128(tmp1, key_ptr[aesKey.nr-1]);
        tmp2 = _mm_aesenc_si128(tmp2, key_ptr[aesKey.nr-1]);
        tmp3 = _mm_aesenc_si128(tmp3, key_ptr[aesKey.nr-1]);
        tmp4 = _mm_aesenc_si128(tmp4, key_ptr[aesKey.nr-1]);
        tmp1 =_mm_aesenclast_si128(tmp1, key_ptr[aesKey.nr]);
        tmp2 =_mm_aesenclast_si128(tmp2, key_ptr[aesKey.nr]);
        tmp3 =_mm_aesenclast_si128(tmp3, key_ptr[aesKey.nr]);
        tmp4 =_mm_aesenclast_si128(tmp4, key_ptr[aesKey.nr]);
        tmp1 = _mm_xor_si128(tmp1, _mm_loadu_si128(&((__m128i*)p.constData())[i*4+0]));
        tmp2 = _mm_xor_si128(tmp2, _mm_loadu_si128(&((__m128i*)p.constData())[i*4+1]));
        tmp3 = _mm_xor_si128(tmp3, _mm_loadu_si128(&((__m128i*)p.constData())[i*4+2]));
        tmp4 = _mm_xor_si128(tmp4, _mm_loadu_si128(&((__m128i*)p.constData())[i*4+3]));

        _mm_storeu_si128(&((__m128i*)c_link)[i*4+0], tmp1);
        _mm_storeu_si128(&((__m128i*)c_link)[i*4+1], tmp2);
        _mm_storeu_si128(&((__m128i*)c_link)[i*4+2], tmp3);
        _mm_storeu_si128(&((__m128i*)c_link)[i*4+3], tmp4);
        tmp1 = _mm_shuffle_epi8(tmp1, BSWAP_MASK);
        tmp2 = _mm_shuffle_epi8(tmp2, BSWAP_MASK);
        tmp3 = _mm_shuffle_epi8(tmp3, BSWAP_MASK);
        tmp4 = _mm_shuffle_epi8(tmp4, BSWAP_MASK);


        X = _mm_xor_si128(X, tmp1);
        X = gfmul_reflected(X, H);
        X = _mm_xor_si128(X, tmp2);
        X = gfmul_reflected(X, H);
        X = _mm_xor_si128(X, tmp3);
        X = gfmul_reflected(X, H);
        X = _mm_xor_si128(X, tmp4);
        X = gfmul_reflected(X, H);
    }
    for(k = i*4; k < p.length()/16; k++){
        tmp1 = _mm_shuffle_epi8(ctr1, BSWAP_EPI64_MASK);
        ctr1 = _mm_add_epi32(ctr1, ONE);
        tmp1 = _mm_xor_si128(tmp1, key_ptr[0]);
        for(j=1; j<aesKey.nr-1; j+=2){
            tmp1 = _mm_aesenc_si128(tmp1, key_ptr[j]);
            tmp1 = _mm_aesenc_si128(tmp1, key_ptr[j+1]);
        }
        tmp1 = _mm_aesenc_si128(tmp1, key_ptr[aesKey.nr-1]);
        tmp1 = _mm_aesenclast_si128(tmp1, key_ptr[aesKey.nr]);
        tmp1 = _mm_xor_si128(tmp1, _mm_loadu_si128(&((__m128i*)p.constData())[k]));
        _mm_storeu_si128(&((__m128i*)c_link)[k], tmp1);
        tmp1 = _mm_shuffle_epi8(tmp1, BSWAP_MASK);
        X =_mm_xor_si128(X, tmp1);
        X = gfmul_reflected(X, H);
    }

    if(p.length() % 16){            // handle last block if necessary
        ctr1 = _mm_shuffle_epi8(ctr1, BSWAP_EPI64_MASK);
        tmp1 = singleAESBlock(ctr1, aesKey);     // encrypt (j0 + 1) -> tmp1

        tmp1 = _mm_xor_si128(tmp1, _mm_loadu_si128(&((__m128i*)p.constData())[k]));
        last_block = tmp1;
        for(j=0; j<p.length() % 16; j++){
            c_link[k*16 + j] = ((unsigned char*) &last_block)[j];
        }
        for(j; j<16; j++){
            ((unsigned char*) &last_block)[j] = 0;
        }
        tmp1 = _mm_shuffle_epi8(last_block, BSWAP_MASK);
        X = _mm_xor_si128(X, tmp1);
        X = gfmul_reflected(X, H);
    }

    // 6. Final Tag
    tmp1 = _mm_insert_epi64(tmp1, p.length()*8, 0);
    tmp1 = _mm_insert_epi64(tmp1, aad.length()*8, 1);

    X = _mm_xor_si128(X, tmp1);
    X = gfmul_reflected(X, H);
    X = _mm_shuffle_epi8(X, BSWAP_MASK);        // bring final GHASH state back to normal world
    T = _mm_xor_si128(X, T);
    _mm_storeu_si128((__m128i*)out.t.data(), T);

    return out;
}


GCM_OUT encrypt_times_four_ghash_times_four(const QByteArray &key, const QByteArray &iv, const QByteArray &aad, const QByteArray &p){
    // 0.1: Check instruction set support
    if(!(AES_GCM_Compatability::hasAES() && AES_GCM_Compatability::hasAVX2() && AES_GCM_Compatability::hasPCLMUL())) {
        throw std::runtime_error("CPU does not support required AES-GCM instructions (AES-NI, PCLMUL, AVX2");
    }

    // 0.2: Check max lengths of paramteres following NIST specification (NIST 800-38d)
    if(key.length() != gcm_keyLength || iv.isEmpty() || iv.length() > MAX_IV_LEN || aad.length() > MAX_AAD_LEN || p.length() > MAX_PLAIN_LEN) {
        return GCM_OUT();
    }

    // 0.3: Key Expansion
    AES_KEY aesKey;
    AES_set_encrypt_key((unsigned char*) key.constData(), gcm_keyLengthBits, &aesKey);
    __m128i* key_ptr = (__m128i*) aesKey.KEY;

    // 0.4 Variable declarations
    __m128i Y0, tmp1, tmp2, tmp3, tmp4, T;
    __m128i ctr1, ctr2, ctr3, ctr4;
    __m128i H, H2, H3, H4;
    __m128i last_block = ZERO;
    __m128i X = ZERO;           // X is used for the running GHASH calculation state
    unsigned int i, j, k;

    GCM_OUT out;
    out.c.resize(p.length());
    out.t.resize(16);
    char* c_link = out.c.data();

    // 1: H = AES_k(0^128) and IV Expansion and T_pre computation
    if(iv.length() == 12){
        Y0 = _mm_loadu_si128((__m128i*) iv.constData());
        Y0 = _mm_insert_epi32(Y0, 0x01000000, 3);       // this is such that in memory a0:| iv0, iv1, iv2, iv3, iv4, ..., iv11, iv12, 00, 00, 00, 01 |:a15

        // encrypt H and Y0
        tmp1 = Y0;
        tmp2 = ZERO;        // future H
        tmp1 = _mm_xor_si128(tmp1, key_ptr[0]);
        tmp2 = _mm_xor_si128(tmp2, key_ptr[0]);
        for(i=1; i < aesKey.nr-1; i+=2){
            tmp1 = _mm_aesenc_si128(tmp1, key_ptr[i]);
            tmp2 = _mm_aesenc_si128(tmp2, key_ptr[i]);
            tmp1 = _mm_aesenc_si128(tmp1, key_ptr[i+1]);
            tmp2 = _mm_aesenc_si128(tmp2, key_ptr[i+1]);
        }
        tmp1 = _mm_aesenc_si128(tmp1, key_ptr[aesKey.nr-1]);
        tmp2 = _mm_aesenc_si128(tmp2, key_ptr[aesKey.nr-1]);
        tmp1 = _mm_aesenclast_si128(tmp1, key_ptr[aesKey.nr]);
        tmp2 = _mm_aesenclast_si128(tmp2, key_ptr[aesKey.nr]);

        H = _mm_shuffle_epi8(tmp2, BSWAP_MASK);
        T = tmp1;

    } else {        // We have to apply GHASH(IV||0^{remainingBytesForFullBlock}||0^32||iv.bitlength())
        // compute H
        tmp2 = _mm_xor_si128(ZERO, key_ptr[0]);
        for(i=1; i<aesKey.nr-1; i+=2){
            tmp2 = _mm_aesenc_si128(tmp2, key_ptr[i]);
            tmp2 = _mm_aesenc_si128(tmp2, key_ptr[i+1]);
        }
        tmp2 = _mm_aesenc_si128(tmp2, key_ptr[aesKey.nr-1]);
        tmp2 = _mm_aesenclast_si128(tmp2, key_ptr[aesKey.nr]);
        H = _mm_shuffle_epi8(tmp2, BSWAP_MASK);

        Y0 = ZERO;
        for(i=0; i < iv.length()/16; i++){  // We do first all full 16 byte blocks
            tmp1 = _mm_loadu_si128(&((__m128i*)iv.constData())[i]);
            tmp1 = _mm_shuffle_epi8(tmp1, BSWAP_MASK);
            Y0 = _mm_xor_si128(Y0, tmp1);
            Y0 = gfmul_reflected(Y0, H);
        }
        if(iv.length() % 16){
            for(j=0; j < iv.length() % 16; j++)
                ((unsigned char*)&last_block)[j] = iv[i*16+j];
            tmp1 = last_block;
            tmp1 = _mm_shuffle_epi8(tmp1, BSWAP_MASK);
            Y0 = _mm_xor_si128(Y0, tmp1);
            Y0 = gfmul_reflected(Y0, H);
        }
        tmp1 = _mm_insert_epi64(tmp1, iv.length()*8, 0);
        tmp1 = _mm_insert_epi64(tmp1, 0, 1);
        Y0 = _mm_xor_si128(Y0, tmp1);
        Y0 = gfmul_reflected(Y0, H);
        Y0 = _mm_shuffle_epi8(Y0, BSWAP_MASK);      // this brings it back into "normal" world

        // compute T_pre
        tmp1 = _mm_xor_si128(Y0, key_ptr[0]);
        for(i=1; i<aesKey.nr-1; i+=2){
            tmp1 = _mm_aesenc_si128(tmp1, key_ptr[i]);
            tmp1 = _mm_aesenc_si128(tmp1, key_ptr[i+1]);
        }
        tmp1 = _mm_aesenc_si128(tmp1, key_ptr[aesKey.nr-1]);
        T = _mm_aesenclast_si128(tmp1, key_ptr[aesKey.nr]);
    }

    // 2. compute constant H keys H2, H3, H4
    H2 = gfmul_reflected(H, H);
    H3 = gfmul_reflected(H2, H);
    H4 = gfmul_reflected(H3, H);

    // 3. GHASH(aad)
    for(i=0; i<aad.length() / 16 / 4; i++){     // apply GHASH on four blocks a time
        tmp1 = _mm_loadu_si128(&((__m128i*)aad.constData())[i*4]);
        tmp2 = _mm_loadu_si128(&((__m128i*)aad.constData())[i*4+1]);
        tmp3 = _mm_loadu_si128(&((__m128i*)aad.constData())[i*4+2]);
        tmp4 = _mm_loadu_si128(&((__m128i*)aad.constData())[i*4+3]);

        tmp1 = _mm_shuffle_epi8(tmp1, BSWAP_MASK);
        tmp2 = _mm_shuffle_epi8(tmp2, BSWAP_MASK);
        tmp3 = _mm_shuffle_epi8(tmp3, BSWAP_MASK);
        tmp4 = _mm_shuffle_epi8(tmp4, BSWAP_MASK);

        tmp1 = _mm_xor_si128(X, tmp1);
        X = gfmul_times_four_reflected(tmp4, tmp3, tmp2, tmp1, H, H2, H3, H4);
    }
    for(i=i*4; i<aad.length()/16; i++){       // apply GHASH on remaining full blocks
        tmp1 = _mm_loadu_si128(&((__m128i*)aad.constData())[i]);
        tmp1 = _mm_shuffle_epi8(tmp1, BSWAP_MASK);
        X = _mm_xor_si128(X, tmp1);
        X = gfmul_reflected(X, H);
    }
    if(aad.length() % 16){                  // apply GHASH on the remaining block if necessary
        last_block = ZERO;
        for(j=0; j<aad.length() % 16; j++){
            ((unsigned char*) &last_block)[j] = aad[i*16 + j];
        }
        tmp1 = _mm_shuffle_epi8(last_block, BSWAP_MASK);
        X = _mm_xor_si128(X, tmp1);
        X = gfmul_reflected(X, H);
    }

    // 4. Ciphertext computation
    ctr1 = _mm_shuffle_epi8(Y0, BSWAP_EPI64_MASK);
    ctr1 = _mm_add_epi32(ctr1, ONE);
    ctr2 = _mm_add_epi32(ctr1, ONE);
    ctr3 = _mm_add_epi32(ctr2, ONE);
    ctr4 = _mm_add_epi32(ctr3, ONE);
    for(i=0; i<p.length() / 16 / 4; i++) {      // process four full blocks at once
        tmp1 = _mm_shuffle_epi8(ctr1, BSWAP_EPI64_MASK);
        tmp2 = _mm_shuffle_epi8(ctr2, BSWAP_EPI64_MASK);
        tmp3 = _mm_shuffle_epi8(ctr3, BSWAP_EPI64_MASK);
        tmp4 = _mm_shuffle_epi8(ctr4, BSWAP_EPI64_MASK);

        ctr1 = _mm_add_epi32(ctr1, FOUR);
        ctr2 = _mm_add_epi32(ctr2, FOUR);
        ctr3 = _mm_add_epi32(ctr3, FOUR);
        ctr4 = _mm_add_epi32(ctr4, FOUR);

        tmp1 =_mm_xor_si128(tmp1, key_ptr[0]);
        tmp2 =_mm_xor_si128(tmp2, key_ptr[0]);
        tmp3 =_mm_xor_si128(tmp3, key_ptr[0]);
        tmp4 =_mm_xor_si128(tmp4, key_ptr[0]);

        for(j=1; j < aesKey.nr-1; j+=2){
            tmp1 = _mm_aesenc_si128(tmp1, key_ptr[j]);
            tmp2 = _mm_aesenc_si128(tmp2, key_ptr[j]);
            tmp3 = _mm_aesenc_si128(tmp3, key_ptr[j]);
            tmp4 = _mm_aesenc_si128(tmp4, key_ptr[j]);

            tmp1 = _mm_aesenc_si128(tmp1, key_ptr[j+1]);
            tmp2 = _mm_aesenc_si128(tmp2, key_ptr[j+1]);
            tmp3 = _mm_aesenc_si128(tmp3, key_ptr[j+1]);
            tmp4 = _mm_aesenc_si128(tmp4, key_ptr[j+1]);
        }

        tmp1 = _mm_aesenc_si128(tmp1, key_ptr[aesKey.nr-1]);
        tmp2 = _mm_aesenc_si128(tmp2, key_ptr[aesKey.nr-1]);
        tmp3 = _mm_aesenc_si128(tmp3, key_ptr[aesKey.nr-1]);
        tmp4 = _mm_aesenc_si128(tmp4, key_ptr[aesKey.nr-1]);

        tmp1 =_mm_aesenclast_si128(tmp1, key_ptr[aesKey.nr]);
        tmp2 =_mm_aesenclast_si128(tmp2, key_ptr[aesKey.nr]);
        tmp3 =_mm_aesenclast_si128(tmp3, key_ptr[aesKey.nr]);
        tmp4 =_mm_aesenclast_si128(tmp4, key_ptr[aesKey.nr]);

        // the _mm_loadu_si128 is necessary because the data might not be 16Byte aligned
        tmp1 = _mm_xor_si128(tmp1, _mm_loadu_si128(&((__m128i*)p.constData())[i*4+0]));
        tmp2 = _mm_xor_si128(tmp2, _mm_loadu_si128(&((__m128i*)p.constData())[i*4+1]));
        tmp3 = _mm_xor_si128(tmp3, _mm_loadu_si128(&((__m128i*)p.constData())[i*4+2]));
        tmp4 = _mm_xor_si128(tmp4, _mm_loadu_si128(&((__m128i*)p.constData())[i*4+3]));

        // the _mm_storeu_si128 is necessary because the destination might not be 16Byte aligned
        _mm_storeu_si128(&((__m128i*)c_link)[i*4+0], tmp1);
        _mm_storeu_si128(&((__m128i*)c_link)[i*4+1], tmp2);
        _mm_storeu_si128(&((__m128i*)c_link)[i*4+2], tmp3);
        _mm_storeu_si128(&((__m128i*)c_link)[i*4+3], tmp4);

        tmp1 = _mm_shuffle_epi8(tmp1, BSWAP_MASK);
        tmp2 = _mm_shuffle_epi8(tmp2, BSWAP_MASK);
        tmp3 = _mm_shuffle_epi8(tmp3, BSWAP_MASK);
        tmp4 = _mm_shuffle_epi8(tmp4, BSWAP_MASK);

        tmp1 = _mm_xor_si128(X, tmp1);
        X = gfmul_times_four_reflected(tmp4, tmp3, tmp2, tmp1, H, H2, H3, H4);
    }

    for(k = i*4; k < p.length()/16; k++){       // handle last full blocks
        tmp1 = _mm_shuffle_epi8(ctr1, BSWAP_EPI64_MASK);
        ctr1 = _mm_add_epi32(ctr1, ONE);
        tmp1 = _mm_xor_si128(tmp1, key_ptr[0]);
        for(j=1; j<aesKey.nr-1; j+=2){
            tmp1 = _mm_aesenc_si128(tmp1, key_ptr[j]);
            tmp1 = _mm_aesenc_si128(tmp1, key_ptr[j+1]);
        }
        tmp1 = _mm_aesenc_si128(tmp1, key_ptr[aesKey.nr-1]);
        tmp1 = _mm_aesenclast_si128(tmp1, key_ptr[aesKey.nr]);
        tmp1 = _mm_xor_si128(tmp1, _mm_loadu_si128(&((__m128i*)p.constData())[k]));
        _mm_storeu_si128(&((__m128i*)c_link)[k], tmp1);
        tmp1 = _mm_shuffle_epi8(tmp1, BSWAP_MASK);
        X =_mm_xor_si128(X, tmp1);
        X = gfmul_reflected(X, H);
    }
    if(p.length() % 16){            // handle last unfull block if necessary
        ctr1 = _mm_shuffle_epi8(ctr1, BSWAP_EPI64_MASK);
        tmp1 = singleAESBlock(ctr1, aesKey);     // encrypt (j0 + 1) -> tmp1

        tmp1 = _mm_xor_si128(tmp1, _mm_loadu_si128(&((__m128i*)p.constData())[k]));
        last_block = tmp1;
        for(j=0; j<p.length() % 16; j++){
            c_link[k*16 + j] = ((unsigned char*) &last_block)[j];
        }
        for(j; j<16; j++){
            ((unsigned char*) &last_block)[j] = 0;
        }
        tmp1 = _mm_shuffle_epi8(last_block, BSWAP_MASK);
        X = _mm_xor_si128(X, tmp1);
        X = gfmul_reflected(X, H);
    }

    // 5. Final Tag
    tmp1 = _mm_insert_epi64(tmp1, p.length()*8, 0);
    tmp1 = _mm_insert_epi64(tmp1, aad.length()*8, 1);

    X = _mm_xor_si128(X, tmp1);
    X = gfmul_reflected(X, H);
    X = _mm_shuffle_epi8(X, BSWAP_MASK);        // bring final GHASH state back to normal world
    T = _mm_xor_si128(X, T);
    _mm_storeu_si128((__m128i*)out.t.data(), T);

    return out;
}

/**
 * @brief encrypt_and_ghash
 * This method is to expected to run inside a thread of which k=4 many should be created.
 * @param threadIdx ranging from [0,3] such that the offsets for data acquisition are clear
 * @param aesKey    pointer to AES_KEY struct
 * @param H_keys    __m128i array containing {H^1, H^2, H^3, H^4}
 * @param Y0
 * @param p         const QByteArray reference. Use p.constData() and load values using _mm_loadu_si128()
 * @param X_ref     passed by value because (even though arithmetically it can not happen), X might change whilst the thread runs this method
 * @param c_link    char* Store values using _mm_storeu_si128()
 * @param X_out     __m128i array, blank, to be filled by X_out[threadIdx] = X
 */
void encrypt_and_ghash_worker(const int32_t& threadIdx, const AES_KEY* const aesKey, const __m128i* const H_keys, const __m128i& Y0, const char* const p_data, const int64_t& p_byte_length, const __m128i X_prev, char* const c_link, __m128i* const X_out){
    __m128i tmp1;
    __m128i X = threadIdx == 0 ? X_prev : ZERO;     // if threadIdx == 0, then start extending X with X_prev as basis
    __m128i* key_ptr = (__m128i*) aesKey->KEY;

    __m128i ctr = _mm_shuffle_epi8(Y0, BSWAP_EPI64_MASK);
    ctr = _mm_add_epi32(ctr, _mm_set_epi32(0, (threadIdx+1), 0, 0));

    for(int i=0; i<p_byte_length / 16 / 4; i++) {
        tmp1 = _mm_shuffle_epi8(ctr, BSWAP_EPI64_MASK);

        ctr = _mm_add_epi32(ctr, FOUR);
        // iv generation is correct

        // encrypt tmp1
        tmp1 =_mm_xor_si128(tmp1, key_ptr[0]);
        for(int j=1; j < aesKey->nr-1; j+=2){
            tmp1 = _mm_aesenc_si128(tmp1, key_ptr[j]);
            tmp1 = _mm_aesenc_si128(tmp1, key_ptr[j+1]);
        }
        tmp1 = _mm_aesenc_si128(tmp1, key_ptr[aesKey->nr-1]);
        tmp1 =_mm_aesenclast_si128(tmp1, key_ptr[aesKey->nr]);

        // the _mm_loadu_si128 is necessary because the data might not be 16Byte aligned
        tmp1 = _mm_xor_si128(tmp1, _mm_loadu_si128(&((__m128i*)p_data)[i*4+threadIdx]));

        // the _mm_storeu_si128 is necessary because the destination might not be 16Byte aligned
        _mm_storeu_si128(&((__m128i*)c_link)[i*4+threadIdx], tmp1);

        // ciphertext generation is also correct

        tmp1 = _mm_shuffle_epi8(tmp1, BSWAP_MASK);
        X = _mm_xor_si128(X, tmp1);

        if(i+1 == (p_byte_length / 16 / 4)){                    // last loop detection
            X = gfmul_reflected(X, H_keys[3-threadIdx]);        //switch(threadIdx){ case 0: H^4; case 1: H^3; case 2: H^2; case 3: H^1 }
        } else{
            X = gfmul_reflected(X, H_keys[3]);                  // else:  H_keys[3] to get H^4
        }
    }

    X_out[threadIdx] = X;
}

GCM_OUT encrypt_times_four_ghash_times_four_parallelized(const QByteArray &key, const QByteArray &iv, const QByteArray &aad, const QByteArray &p){
    // 0.1: Check instruction set support
    if(!(AES_GCM_Compatability::hasAES() && AES_GCM_Compatability::hasAVX2() && AES_GCM_Compatability::hasPCLMUL())) {
        throw std::runtime_error("CPU does not support required AES-GCM instructions (AES-NI, PCLMUL, AVX2");
    }

    // 0.2: Check max lengths of paramteres following NIST specification (NIST 800-38d)
    if(key.length() != gcm_keyLength || iv.isEmpty() || iv.length() > MAX_IV_LEN || aad.length() > MAX_AAD_LEN || p.length() > MAX_PLAIN_LEN) {
        return GCM_OUT();
    }

    // 0.3: Key Expansion
    AES_KEY aesKey;
    AES_set_encrypt_key((unsigned char*) key.constData(), gcm_keyLengthBits, &aesKey);
    __m128i* key_ptr = (__m128i*) aesKey.KEY;

    // 0.4 Variable declarations
    __m128i Y0, tmp1, tmp2, tmp3, tmp4, T;
    __m128i ctr1;
    __m128i H, H2, H3, H4;
    __m128i last_block = ZERO;
    __m128i X = ZERO;           // X is used for the running GHASH calculation state
    unsigned int i, j, k;

    GCM_OUT out;
    out.c.resize(p.length());
    out.t.resize(16);
    char* c_link = out.c.data();

    // 1: H = AES_k(0^128) and IV Expansion and T_pre computation
    if(iv.length() == 12){
        Y0 = _mm_loadu_si128((__m128i*) iv.constData());
        Y0 = _mm_insert_epi32(Y0, 0x01000000, 3);       // this is such that in memory a0:| iv0, iv1, iv2, iv3, iv4, ..., iv11, iv12, 00, 00, 00, 01 |:a15

        // encrypt H and Y0
        tmp1 = Y0;
        tmp2 = ZERO;        // future H
        tmp1 = _mm_xor_si128(tmp1, key_ptr[0]);
        tmp2 = _mm_xor_si128(tmp2, key_ptr[0]);
        for(i=1; i < aesKey.nr-1; i+=2){
            tmp1 = _mm_aesenc_si128(tmp1, key_ptr[i]);
            tmp2 = _mm_aesenc_si128(tmp2, key_ptr[i]);
            tmp1 = _mm_aesenc_si128(tmp1, key_ptr[i+1]);
            tmp2 = _mm_aesenc_si128(tmp2, key_ptr[i+1]);
        }
        tmp1 = _mm_aesenc_si128(tmp1, key_ptr[aesKey.nr-1]);
        tmp2 = _mm_aesenc_si128(tmp2, key_ptr[aesKey.nr-1]);
        tmp1 = _mm_aesenclast_si128(tmp1, key_ptr[aesKey.nr]);
        tmp2 = _mm_aesenclast_si128(tmp2, key_ptr[aesKey.nr]);

        H = _mm_shuffle_epi8(tmp2, BSWAP_MASK);
        T = tmp1;

    } else {        // We have to apply GHASH(IV||0^{remainingBytesForFullBlock}||0^32||iv.bitlength())
        // compute H
        tmp2 = _mm_xor_si128(ZERO, key_ptr[0]);
        for(i=1; i<aesKey.nr-1; i+=2){
            tmp2 = _mm_aesenc_si128(tmp2, key_ptr[i]);
            tmp2 = _mm_aesenc_si128(tmp2, key_ptr[i+1]);
        }
        tmp2 = _mm_aesenc_si128(tmp2, key_ptr[aesKey.nr-1]);
        tmp2 = _mm_aesenclast_si128(tmp2, key_ptr[aesKey.nr]);
        H = _mm_shuffle_epi8(tmp2, BSWAP_MASK);

        Y0 = ZERO;
        for(i=0; i < iv.length()/16; i++){  // We do first all full 16 byte blocks
            tmp1 = _mm_loadu_si128(&((__m128i*)iv.constData())[i]); // the _mm_loadu_si128 is necessary because the data might not be 16Byte aligned
            tmp1 = _mm_shuffle_epi8(tmp1, BSWAP_MASK);
            Y0 = _mm_xor_si128(Y0, tmp1);
            Y0 = gfmul_reflected(Y0, H);
        }
        if(iv.length() % 16){
            for(j=0; j < iv.length() % 16; j++)
                ((unsigned char*)&last_block)[j] = iv[i*16+j];
            tmp1 = last_block;
            tmp1 = _mm_shuffle_epi8(tmp1, BSWAP_MASK);
            Y0 = _mm_xor_si128(Y0, tmp1);
            Y0 = gfmul_reflected(Y0, H);
        }
        tmp1 = _mm_insert_epi64(tmp1, iv.length()*8, 0);
        tmp1 = _mm_insert_epi64(tmp1, 0, 1);
        Y0 = _mm_xor_si128(Y0, tmp1);
        Y0 = gfmul_reflected(Y0, H);
        Y0 = _mm_shuffle_epi8(Y0, BSWAP_MASK);      // this brings it back into "normal" world

        // compute T_pre
        tmp1 = _mm_xor_si128(Y0, key_ptr[0]);
        for(i=1; i<aesKey.nr-1; i+=2){
            tmp1 = _mm_aesenc_si128(tmp1, key_ptr[i]);
            tmp1 = _mm_aesenc_si128(tmp1, key_ptr[i+1]);
        }
        tmp1 = _mm_aesenc_si128(tmp1, key_ptr[aesKey.nr-1]);
        T = _mm_aesenclast_si128(tmp1, key_ptr[aesKey.nr]);
    }

    // 2. compute constant H keys H2, H3, H4
    H2 = gfmul_reflected(H, H);
    H3 = gfmul_reflected(H2, H);
    H4 = gfmul_reflected(H3, H);

    // 3. GHASH(aad) (for starters, we do not parallelize it, as aad is normally not so big. One for-loop handles 64Bytes of aad, only if aad is >>> k*64Bytes, then a parallelization would be handy)
    for(i=0; i<aad.length() / 16 / 4; i++){     // apply GHASH on four blocks a time
        // the _mm_loadu_si128 is necessary because the data might not be 16Byte aligned
        tmp1 = _mm_loadu_si128(&((__m128i*)aad.constData())[i*4]);
        tmp2 = _mm_loadu_si128(&((__m128i*)aad.constData())[i*4+1]);
        tmp3 = _mm_loadu_si128(&((__m128i*)aad.constData())[i*4+2]);
        tmp4 = _mm_loadu_si128(&((__m128i*)aad.constData())[i*4+3]);

        tmp1 = _mm_shuffle_epi8(tmp1, BSWAP_MASK);
        tmp2 = _mm_shuffle_epi8(tmp2, BSWAP_MASK);
        tmp3 = _mm_shuffle_epi8(tmp3, BSWAP_MASK);
        tmp4 = _mm_shuffle_epi8(tmp4, BSWAP_MASK);

        tmp1 = _mm_xor_si128(X, tmp1);
        X = gfmul_times_four_reflected(tmp4, tmp3, tmp2, tmp1, H, H2, H3, H4);
    }
    for(i=i*4; i<aad.length()/16; i++){       // apply GHASH on remaining full blocks
        tmp1 = _mm_loadu_si128(&((__m128i*)aad.constData())[i]);
        tmp1 = _mm_shuffle_epi8(tmp1, BSWAP_MASK);
        X = _mm_xor_si128(X, tmp1);
        X = gfmul_reflected(X, H);
    }
    if(aad.length() % 16){                  // apply GHASH on the remaining block if necessary
        last_block = ZERO;
        for(j=0; j<aad.length() % 16; j++){
            ((unsigned char*) &last_block)[j] = aad[i*16 + j];
        }
        tmp1 = _mm_shuffle_epi8(last_block, BSWAP_MASK);
        X = _mm_xor_si128(X, tmp1);
        X = gfmul_reflected(X, H);
    }

    // 4. Ciphertext computation
    // ------------
    // Parallelization:
    // The idea, we use the Parallel GHASH routine as defined in Intel Doc. B Figure 8
    // While their illustration only covers the cases for n % 4 == 0, if we have more blocks we just add them normally (sequentially)
    // this parallelization is fixed to 4 simultanous running threads. The choice for k=4, should maybe later be dynimacally made dependent on the number of cores available with AVX2 instruction set.


    /*
    for(int m=0; m<4; m++){
        ctr1 = _mm_shuffle_epi8(Y0, BSWAP_EPI64_MASK);
        ctr1 = _mm_add_epi32(ctr1, ONE);
        ctr2 = _mm_add_epi32(ctr1, ONE);
        ctr3 = _mm_add_epi32(ctr2, ONE);
        ctr4 = _mm_add_epi32(ctr3, ONE);

        for(i=0; i<p.length() / 16 / 4; i++) {      // process four full blocks at once
            qInfo()<<"in loop";
            tmp1 = _mm_shuffle_epi8(ctr1, BSWAP_EPI64_MASK);
            tmp2 = _mm_shuffle_epi8(ctr2, BSWAP_EPI64_MASK);
            tmp3 = _mm_shuffle_epi8(ctr3, BSWAP_EPI64_MASK);
            tmp4 = _mm_shuffle_epi8(ctr4, BSWAP_EPI64_MASK);
            qInfo()<<"[encrypt_times_four_ghash_times_four_parallelized]"<<"ctr1"<<": "<<print128_hex_lanes(tmp1);
            qInfo()<<"[encrypt_times_four_ghash_times_four_parallelized]"<<"ctr2"<<": "<<print128_hex_lanes(tmp2);
            qInfo()<<"[encrypt_times_four_ghash_times_four_parallelized]"<<"ctr3"<<": "<<print128_hex_lanes(tmp3);
            qInfo()<<"[encrypt_times_four_ghash_times_four_parallelized]"<<"ctr4"<<": "<<print128_hex_lanes(tmp4);

            ctr1 = _mm_add_epi32(ctr1, FOUR);
            ctr2 = _mm_add_epi32(ctr2, FOUR);
            ctr3 = _mm_add_epi32(ctr3, FOUR);
            ctr4 = _mm_add_epi32(ctr4, FOUR);

            tmp1 =_mm_xor_si128(tmp1, key_ptr[0]);
            tmp2 =_mm_xor_si128(tmp2, key_ptr[0]);
            tmp3 =_mm_xor_si128(tmp3, key_ptr[0]);
            tmp4 =_mm_xor_si128(tmp4, key_ptr[0]);

            for(j=1; j < aesKey.nr-1; j+=2){
                tmp1 = _mm_aesenc_si128(tmp1, key_ptr[j]);
                tmp2 = _mm_aesenc_si128(tmp2, key_ptr[j]);
                tmp3 = _mm_aesenc_si128(tmp3, key_ptr[j]);
                tmp4 = _mm_aesenc_si128(tmp4, key_ptr[j]);

                tmp1 = _mm_aesenc_si128(tmp1, key_ptr[j+1]);
                tmp2 = _mm_aesenc_si128(tmp2, key_ptr[j+1]);
                tmp3 = _mm_aesenc_si128(tmp3, key_ptr[j+1]);
                tmp4 = _mm_aesenc_si128(tmp4, key_ptr[j+1]);
            }

            tmp1 = _mm_aesenc_si128(tmp1, key_ptr[aesKey.nr-1]);
            tmp2 = _mm_aesenc_si128(tmp2, key_ptr[aesKey.nr-1]);
            tmp3 = _mm_aesenc_si128(tmp3, key_ptr[aesKey.nr-1]);
            tmp4 = _mm_aesenc_si128(tmp4, key_ptr[aesKey.nr-1]);

            tmp1 =_mm_aesenclast_si128(tmp1, key_ptr[aesKey.nr]);
            tmp2 =_mm_aesenclast_si128(tmp2, key_ptr[aesKey.nr]);
            tmp3 =_mm_aesenclast_si128(tmp3, key_ptr[aesKey.nr]);
            tmp4 =_mm_aesenclast_si128(tmp4, key_ptr[aesKey.nr]);

            __m128i tt[4] = {tmp1, tmp2, tmp3, tmp4};

            // the _mm_loadu_si128 is necessary because the data might not be 16Byte aligned
            tt[m] = _mm_xor_si128(tt[m], _mm_loadu_si128(&((__m128i*)p.constData())[i*4+m]));
            //tmp2 = _mm_xor_si128(tmp2, _mm_loadu_si128(&((__m128i*)p.constData())[i*4+1]));
            //tmp3 = _mm_xor_si128(tmp3, _mm_loadu_si128(&((__m128i*)p.constData())[i*4+2]));
            //tmp4 = _mm_xor_si128(tmp4, _mm_loadu_si128(&((__m128i*)p.constData())[i*4+3]));
            qInfo()<<"[encrypt_times_four_ghash_times_four_parallelized]"<<"c"<<i<<"- "<<m<<":"<<print128_hex_lanes(tt[m]);
            //qInfo()<<"[encrypt_times_four_ghash_times_four_parallelized]"<<"c"<<i<<"- 1:"<<print128_hex_lanes(tmp2);
            //qInfo()<<"[encrypt_times_four_ghash_times_four_parallelized]"<<"c"<<i<<"- 2:"<<print128_hex_lanes(tmp3);
            //qInfo()<<"[encrypt_times_four_ghash_times_four_parallelized]"<<"c"<<i<<"- 0:"<<print128_hex_lanes(tmp4);

            // the _mm_storeu_si128 is necessary because the destination might not be 16Byte aligned
            _mm_storeu_si128(&((__m128i*)c_link)[i*4+m], tt[m]);
            //_mm_storeu_si128(&((__m128i*)c_link)[i*4+1], tmp2);
            //_mm_storeu_si128(&((__m128i*)c_link)[i*4+2], tmp3);
            //_mm_storeu_si128(&((__m128i*)c_link)[i*4+3], tmp4);

            tmp1 = _mm_shuffle_epi8(tmp1, BSWAP_MASK);
            tmp2 = _mm_shuffle_epi8(tmp2, BSWAP_MASK);
            tmp3 = _mm_shuffle_epi8(tmp3, BSWAP_MASK);
            tmp4 = _mm_shuffle_epi8(tmp4, BSWAP_MASK);

            tmp1 = _mm_xor_si128(X, tmp1);
            X = gfmul_times_four_reflected(tmp4, tmp3, tmp2, tmp1, H, H2, H3, H4);
        }
    }
    */



    int num_threads = 4;        // std::thread::hardware_capability
    __m128i H_keys[4] = {H, H2, H3, H4};   // must be adapted to num_threads
    __m128i X_out[num_threads];
    //std::thread threads[num_threads];
    qInfo()<<"P byte length"<<p.length();
    for(int m=0; m<num_threads; m++){
        encrypt_and_ghash_worker(m, &aesKey, (__m128i*) H_keys, Y0, p.constData(), p.length(), X, c_link, (__m128i*) X_out);
        //threads[m] = std::thread(encrypt_and_ghash_worker, m, &aesKey, H_keys, Y0, p.constData(), p.length(), X, c_link, (__m128i*) X_out);
    }

    for(int m=0; m<num_threads; m++){
        //threads[m].join();
        X = _mm_xor_si128(X, X_out[m]);
    }



    // The code below is the same thing:
    //ctr1 = _mm_add_epi32(ctr1, ONE);
    //for(i=0; i<p.length() / 16 / 4; i++) ctr1 = _mm_add_epi32(ctr1, FOUR);
    i = p.length()/16/4;
    ctr1 = _mm_shuffle_epi8(Y0, BSWAP_EPI64_MASK);
    ctr1 = _mm_add_epi32(ctr1, _mm_set_epi32(0, 1 + 4*i, 0, 0));

    for(k = i*4; k < p.length()/16; k++){       // handle last full blocks
        tmp1 = _mm_shuffle_epi8(ctr1, BSWAP_EPI64_MASK);
        ctr1 = _mm_add_epi32(ctr1, ONE);
        tmp1 = _mm_xor_si128(tmp1, key_ptr[0]);
        for(j=1; j<aesKey.nr-1; j+=2){
            tmp1 = _mm_aesenc_si128(tmp1, key_ptr[j]);
            tmp1 = _mm_aesenc_si128(tmp1, key_ptr[j+1]);
        }
        tmp1 = _mm_aesenc_si128(tmp1, key_ptr[aesKey.nr-1]);
        tmp1 = _mm_aesenclast_si128(tmp1, key_ptr[aesKey.nr]);
        tmp1 = _mm_xor_si128(tmp1, _mm_loadu_si128(&((__m128i*)p.constData())[k]));
        _mm_storeu_si128(&((__m128i*)c_link)[k], tmp1); // the _mm_storeu_si128 is necessary because the destination might not be 16Byte aligned
        tmp1 = _mm_shuffle_epi8(tmp1, BSWAP_MASK);
        X =_mm_xor_si128(X, tmp1);
        X = gfmul_reflected(X, H);
    }
    if(p.length() % 16){            // handle last unfull block if necessary
        ctr1 = _mm_shuffle_epi8(ctr1, BSWAP_EPI64_MASK);
        tmp1 = singleAESBlock(ctr1, aesKey);     // encrypt (j0 + 1) -> tmp1

        tmp1 = _mm_xor_si128(tmp1, _mm_loadu_si128(&((__m128i*)p.constData())[k]));
        last_block = tmp1;
        for(j=0; j<p.length() % 16; j++){
            c_link[k*16 + j] = ((unsigned char*) &last_block)[j];
        }
        for(j; j<16; j++){
            ((unsigned char*) &last_block)[j] = 0;
        }
        tmp1 = _mm_shuffle_epi8(last_block, BSWAP_MASK);
        X = _mm_xor_si128(X, tmp1);
        X = gfmul_reflected(X, H);
    }

    // 5. Final Tag
    tmp1 = _mm_insert_epi64(tmp1, p.length()*8, 0);
    tmp1 = _mm_insert_epi64(tmp1, aad.length()*8, 1);

    X = _mm_xor_si128(X, tmp1);
    X = gfmul_reflected(X, H);
    X = _mm_shuffle_epi8(X, BSWAP_MASK);        // bring final GHASH state back to normal world
    T = _mm_xor_si128(X, T);
    _mm_storeu_si128((__m128i*)out.t.data(), T);

    return out;
}



#include <QString>
/** FOR ASSERTION VALUES
 * You must load them in such that for 'feffe9928665731c6d6a8f9467308308' corresponds to {fe, ff, e9, ..., 83, 08}
 * You can not use _mm_set_epi32, as that will reverse the order
 * @brief gcm_test
 */

int add(int i){
    return i+1;
}

void gcm_test(){
    // This is Test 4 from revised NIST
    QByteArray key = QByteArray::fromHex("feffe9928665731c6d6a8f9467308308");       // below is the same operation
    //unsigned char key[16] = {0xfe, 0xff, 0xe9, 0x92, 0x86, 0x65, 0x73, 0x1c, 0x6d, 0x6a, 0x8f, 0x94, 0x67, 0x30, 0x83, 0x08};
    QByteArray p =  QByteArray::fromHex("d9313225f88406e5a55909c5aff5269a86a7a9531534f7da2e4c303d8a318a721c3c0c95956809532fcf0e2449a6b525b16aedf5aa0de657ba637b39");
    QByteArray iv = QByteArray::fromHex("cafebabefacedbaddecaf888");        // will result in {ca, fe, ba, be, ..., f8, 88} array, with arr[0] = ca
    QByteArray a = QByteArray::fromHex("feedfacedeadbeeffeedfacedeadbeefabaddad2");

    QByteArray c_assert = QByteArray::fromHex("42831ec2217774244b7221b784d0d49ce3aa212f2c02a4e035c17e2329aca12e21d514b25466931c7d8f6a5aac84aa051ba30b396a0aac973d58e091");
    QByteArray t_assert = QByteArray::fromHex("5bc94fbc3221a5db94fae95ae7121a47");

    GCM_OUT out_test = encrypt(key, iv, a, p);
    if(c_assert == out_test.c){
        qInfo() << "Assertion ( \"C\" ): OK";
    } else{
        qInfo() << "Assertion ( \"C\" ): WRONG";
    }
    if(t_assert == out_test.t){
        qInfo() << "Assertion ( \"T\" ): OK";
    } else{
        qInfo() << "Assertion ( \"T\" ): WRONG";
    }

    GCM_OUT out_test2 = encrypt_times_four(key, iv, a, p);
    if(c_assert == out_test2.c){
        qInfo() << "Assertion ( \"C\" ): OK";
    } else{
        qInfo() << "Assertion ( \"C\" ): WRONG";
    }
    if(t_assert == out_test2.t){
        qInfo() << "Assertion ( \"T\" ): OK";
    } else{
        qInfo() << "Assertion ( \"T\" ): WRONG";
    }

    GCM_OUT out_test3 = encrypt_times_four_ghash_times_four(key, iv, a, p);
    if(c_assert == out_test3.c){
        qInfo() << "Assertion ( \"C\" ): OK";
    } else{
        qInfo() << "Assertion ( \"C\" ): WRONG";
    }
    if(t_assert == out_test3.t){
        qInfo() << "Assertion ( \"T\" ): OK";
    } else{
        qInfo() << "Assertion ( \"T\" ): WRONG";
    }

    GCM_OUT out_test4 = encrypt_times_four_ghash_times_four_parallelized(key, iv, a, p);
    if(c_assert == out_test4.c){
        qInfo() << "Assertion ( \"C\" ): OK";
    } else{
        qInfo() << "Assertion ( \"C\" ): WRONG";
    }
    if(t_assert == out_test4.t){
        qInfo() << "Assertion ( \"T\" ): OK";
    } else{
        qInfo() << "Assertion ( \"T\" ): WRONG";
    }
    qInfo()<<"---------------------------";

    // This is Test 6 from revised NIST
    key = QByteArray::fromHex("feffe9928665731c6d6a8f9467308308");
    p = QByteArray::fromHex("d9313225f88406e5a55909c5aff5269a86a7a9531534f7da2e4c303d8a318a721c3c0c95956809532fcf0e2449a6b525b16aedf5aa0de657ba637b39");
    a = QByteArray::fromHex("feedfacedeadbeeffeedfacedeadbeefabaddad2");
    iv = QByteArray::fromHex("9313225df88406e555909c5aff5269aa6a7a9538534f7da1e4c303d2a318a728c3c0c95156809539fcf0e2429a6b525416aedbf5a0de6a57a637b39b");
    c_assert = QByteArray::fromHex("8ce24998625615b603a033aca13fb894be9112a5c3a211a8ba262a3cca7e2ca701e4a9a4fba43c90ccdcb281d48c7c6fd62875d2aca417034c34aee5");
    t_assert = QByteArray::fromHex("619cc5aefffe0bfa462af43c1699d050");

    out_test = encrypt(key, iv, a, p);
    if(c_assert == out_test.c){
        qInfo() << "Assertion ( \"C\" ): OK";
    } else{
        qInfo() << "Assertion ( \"C\" ): WRONG";
    }
    if(t_assert == out_test.t){
        qInfo() << "Assertion ( \"T\" ): OK";
    } else{
        qInfo() << "Assertion ( \"T\" ): WRONG";
    }

    out_test2 = encrypt_times_four(key, iv, a, p);
    if(c_assert == out_test2.c){
        qInfo() << "Assertion ( \"C\" ): OK";
    } else{
        qInfo() << "Assertion ( \"C\" ): WRONG";
    }
    if(t_assert == out_test2.t){
        qInfo() << "Assertion ( \"T\" ): OK";
    } else{
        qInfo() << "Assertion ( \"T\" ): WRONG";
    }

    out_test3 = encrypt_times_four_ghash_times_four(key, iv, a, p);
    if(c_assert == out_test3.c){
        qInfo() << "Assertion ( \"C\" ): OK";
    } else{
        qInfo() << "Assertion ( \"C\" ): WRONG";
    }
    if(t_assert == out_test3.t){
        qInfo() << "Assertion ( \"T\" ): OK";
    } else{
        qInfo() << "Assertion ( \"T\" ): WRONG";
    }

    out_test4 = encrypt_times_four_ghash_times_four_parallelized(key, iv, a, p);
    if(c_assert == out_test4.c){
        qInfo() << "Assertion ( \"C\" ): OK";
    } else{
        qInfo() << "Assertion ( \"C\" ): WRONG";
    }
    if(t_assert == out_test4.t){
        qInfo() << "Assertion ( \"T\" ): OK";
    } else{
        qInfo() << "Assertion ( \"T\" ): WRONG";
    }
    qInfo()<<"---------------------------";


    // This is custom test
    key = QByteArray::fromHex("58ffe9928ab5731c6d6a8f9467308314");
    p = QByteArray::fromHex("d9313225f88406e5a55909c5aff5269a86a7a9531534f7da2e4c303d8a318a721c3c0c95956809532fcf0e2449a6b525b16aedf5aa0de657ba637b39a8b4690145d3feea153e840bc80c95956809532fcf0e2449a6b525b16aedf5aa0de657ba637b39a8b4690145d3feea1d9313225f88406e1847abc630ad1388efdb9d0aabc9e01305063224ebfdbd810c7a0d041431405a55909c5aff5269a86a7a95837cd83a0341df9e93132e4bf90ca7d9e921b9401ffe8abaaccade310981dace325f88406e5a55909c5aff5269a86a7ad8adebe83988444be7fa818460c0edbe87839014123a9531534f7da2e4c303d8a318a721c3c0c95956809532fcf0e2449a6b525b169531534f7da2e4c303d8a311e4c303d2a318a7");
    a = QByteArray::fromHex("58e32f3b790b0983ed1da18ef020f6bc2edb790acdcf31be8ad573667d678557edbdd56b94c56ee111fb3fb862ffc391b90dc0c4d95656fbd545bb8cf1833aa315665963f99863b4c01adde03b26038f10112c1cc6a881d21c6d09479ecf75dd76b739c806328cc77cddf304b65fddd5bba51c21a0887c9e59a4fa1519672ad205c9e50f3d5f0536bf3b60b7dd3ba4f7789d9aaf94812b095151d0");
    iv = QByteArray::fromHex("9313225df88406e555909c5aff5269aa6a7a9538534f7da1e4c303d2a318a728c3c0c95156809539fcf0e2429a6b525416aedbf5a0de6a57a637b39b");
    c_assert = QByteArray::fromHex("b397e284168442e5425c8437a170ec442916a87529c3507897fd384e76f0902c0e55ef073a958ad204a2da4fe8545217c87a7c1010d7c0cdca670ab53638dd8f9e3ba1fae9e528dcefa8010d5d23ec58c0c827458774ad2516cc1642b51c9dba6b0c297f4bb524dade1241ff050c7fe2ed7fe7bbde323539dab4ec695c3818ea01a4aa52a35c9434fc936223f6f790adc70917a5df70b1113150d41c3290a9474f48e5993b0464f350cc81d49a48e70b2e3fc2302bae6b58957e952aecdb942c51585809824d9ebb67eb24362f0472c080b18eaefb27e786e8ec05d2168bc163a66c6d582ea234eff819f5f712b93b3f798dae3d6b0c173c135f1463424bca4d6a9100e14d087753c750f6ebe318f4bb60d7011104c9");
    t_assert = QByteArray::fromHex("921632a41d9faf02cf8bb48a8fc150de");

    out_test = encrypt(key, iv, a, p);
    if(c_assert == out_test.c){
        qInfo() << "Assertion ( \"C\" ): OK";
    } else{
        qInfo() << "Assertion ( \"C\" ): WRONG";
    }
    if(t_assert == out_test.t){
        qInfo() << "Assertion ( \"T\" ): OK";
    } else{
        qInfo() << "Assertion ( \"T\" ): WRONG";
    }

    out_test2 = encrypt_times_four(key, iv, a, p);
    if(c_assert == out_test2.c){
        qInfo() << "Assertion ( \"C\" ): OK";
    } else{
        qInfo() << "Assertion ( \"C\" ): WRONG";
    }
    if(t_assert == out_test2.t){
        qInfo() << "Assertion ( \"T\" ): OK";
    } else{
        qInfo() << "Assertion ( \"T\" ): WRONG";
    }

    out_test3 = encrypt_times_four_ghash_times_four(key, iv, a, p);
    if(c_assert == out_test3.c){
        qInfo() << "Assertion ( \"C\" ): OK";
    } else{
        qInfo() << "Assertion ( \"C\" ): WRONG";
    }
    if(t_assert == out_test3.t){
        qInfo() << "Assertion ( \"T\" ): OK";
    } else{
        qInfo() << "Assertion ( \"T\" ): WRONG";
    }

    out_test4 = encrypt_times_four_ghash_times_four_parallelized(key, iv, a, p);
    if(c_assert == out_test4.c){
        qInfo() << "Assertion ( \"C\" ): OK";
    } else{
        qInfo() << "Assertion ( \"C\" ): WRONG";
    }
    if(t_assert == out_test4.t){
        qInfo() << "Assertion ( \"T\" ): OK";
    } else{
        qInfo() << "Assertion ( \"T\" ): WRONG";
    }
    qInfo()<<"---------------------------";

    QByteArray iv_t = BenchmarkUtil::randQByteArray(128, 96);
    QByteArray a_t = BenchmarkUtil::randQByteArray(500, 250);
    QByteArray p_t = BenchmarkUtil::randQByteArray(2000, 700);
    //BenchmarkUtil::runEncryptBenchmark("Normal_Encrypt_Single", encrypt, key, iv_t, a_t, p_t);
    //BenchmarkUtil::runEncryptBenchmark("Normal_Encrypt_Times_Four", encrypt_times_four, key, iv_t, a_t, p_t);
    //BenchmarkUtil::runEncryptBenchmark("Reduce_Four_Encrypt_Times_Four", encrypt_times_four_ghash_times_four, key, iv_t, a_t, p_t);
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
