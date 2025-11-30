#include <gcm/gcm.h>
#include <gcm/gfmul.h>
#include <gcm/gcm_compatability.h>
#include <aesni/aesni-key-init.h>

#include <smmintrin.h>
#include <qaesencryption.h>
#include <QDebug>
#include <QString>

void singleAESBlock(const unsigned char* in, unsigned char* out, const unsigned char* key, int number_of_rounds){
    __m128i tmp = _mm_loadu_si128((__m128i*) in);   // take first 16 bytes from in and store in tmp pointer
    tmp = _mm_xor_si128(tmp, ((__m128i*)key)[0]);   // apply Round 0 key
    for(int i=1; i<number_of_rounds; i++){
        tmp = _mm_aesenc_si128(tmp, ((__m128i*)key)[i]);    // Rounds: [1 - (n-1)]
    }
    tmp = _mm_aesenclast_si128(tmp, ((__m128i*)key)[number_of_rounds]); // Round: n
    _mm_storeu_si128((__m128i*)out, tmp);           // store tmp in first 16 bytes of out pointer
}

__m128i singleAESBlock(const __m128i& in, const __m128i* const key, int number_of_rounds){
    __m128i out = in;
    out = _mm_xor_si128(out, key[0]);
    for(int i=1; i<number_of_rounds; i++){
        out = _mm_aesenc_si128(out, key[i]);
    }
    out = _mm_aesenclast_si128(out, key[number_of_rounds]);
    return out;
}

struct GCM_OUT encrypt(const QByteArray &key, const QByteArray &iv, const QByteArray &aad, const QByteArray &p){
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

    // 1: H = AES_k(0^128)
    __m128i H = singleAESBlock(ZERO, (__m128i*) aesKey.KEY, gcm_n_rounds);

    // 2. IV Expansion
    __m128i j0;
    if(iv.length() == 12){
        j0 = _mm_loadu_si128((__m128i*) iv.constData());
        j0 = _mm_insert_epi32(j0, 0x01000000, 3);       // this is such that in memory a0:| iv0, iv1, iv2, iv3, iv4, ..., iv11, iv12, 00, 00, 00, 01 |:a15
    } else {

    }



    return GCM_OUT();

}

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


#include <QString>
/** FOR ASSERTION VALUES
 * You must load them in such that for 'feffe9928665731c6d6a8f9467308308' corresponds to {fe, ff, e9, ..., 83, 08}
 * You can not use _mm_set_epi32, as that will reverse the order
 * @brief gcm_test
 */


void gcm_test(){
    //feffe9928665731c6d6a8f9467308308
    QByteArray key = QByteArray::fromHex("feffe9928665731c6d6a8f9467308308");       // below is the same operation
    //unsigned char key[16] = {0xfe, 0xff, 0xe9, 0x92, 0x86, 0x65, 0x73, 0x1c, 0x6d, 0x6a, 0x8f, 0x94, 0x67, 0x30, 0x83, 0x08};
    QByteArray p =  QByteArray::fromHex("d9313225f88406e5a55909c5aff5269a86a7a9531534f7da2e4c303d8a318a721c3c0c95956809532fcf0e2449a6b525b16aedf5aa0de657ba637b39");
    QByteArray iv = QByteArray::fromHex("cafebabefacedbaddecaf888");        // will result in {ca, fe, ba, be, ..., f8, 88} array, with arr[0] = ca
    //qInfo()<< "iv a11:|"<<QString::asprintf("%02X %02X %02X %02X %02X %02X %02X %02X %02X %02X %02X %02X", iv[11], iv[10], iv[9], iv[8], iv[7], iv[6], iv[5], iv[4], iv[3], iv[2], iv[1], iv[0])<<"|:a0";
    QByteArray a = QByteArray::fromHex("feedfacedeadbeeffeedfacedeadbeefabaddad2");
    unsigned char out[p.length()];
    unsigned char tag[16];


    // 0.3: Key Expansion
    AES_KEY aesKey;
    AES_set_encrypt_key((unsigned char*) key.constData(), gcm_keyLengthBits, &aesKey);
    __m128i KEY[gcm_n_rounds+1];
    memcpy((unsigned char*) KEY, aesKey.KEY, 16*(gcm_n_rounds+1));

    // 1: H = AES_k(0^128)
    __m128i H = singleAESBlock(ZERO, KEY, gcm_n_rounds);
    H = _mm_shuffle_epi8(H, BSWAP_MASK);
    //qInfo()<<"H: a15:|"<<print128_hex_lanes(H)<<"|:a0";
    __m128i h_assert = _mm_loadu_si128((__m128i*) QByteArray::fromHex("b83b533708bf535d0aa6e52980d53b78").constData());
    assert_m128(_mm_shuffle_epi8(H, BSWAP_MASK), h_assert, "H");



    // 2. IV Expansion to J0
    __m128i j0;
    if(iv.length() == 12){
        j0 = _mm_loadu_si128((__m128i*) iv.constData());
        j0 = _mm_insert_epi32(j0, 0x01000000, 3);       // this is such that in memory a0:| iv0, iv1, iv2, iv3, iv4, ..., iv11, iv12, 00, 00, 00, 01 |:a15
    } else {
        qInfo()<<"iv not length matching for test";
    }
    qInfo()<<"J0 a15:|"<<print128_hex_lanes(j0)<<"|:a0";
    __m128i j0_assert = _mm_loadu_si128((__m128i*) QByteArray::fromHex("cafebabefacedbaddecaf88800000001").constData());
    assert_m128(j0, j0_assert, "J0");

    // 2.1 Tag_key = AES_k(J0)
    // create Tag T
    __m128i T = singleAESBlock(j0, KEY, gcm_n_rounds);
    __m128i T_pre_assert = _mm_loadu_si128((__m128i*)QByteArray::fromHex("3247184b3c4f69a44dbcd22887bbb418").constData());
    assert_m128(T, T_pre_assert, "T_pre");



    int i,j;
    __m128i tmp1, tmp2, tmp3, tmp4, ctr1, ctr2, ctr3, ctr4;
    __m128i X = _mm_setzero_si128();
    __m128i X_asserts[6] = {_mm_loadu_si128((__m128i*) QByteArray::fromHex("ed56aaf8a72d67049fdb9228edba1322").constData()),
                            _mm_loadu_si128((__m128i*) QByteArray::fromHex("cd47221ccef0554ee4bb044c88150352").constData()),
                            _mm_loadu_si128((__m128i*) QByteArray::fromHex("54f5e1b2b5a8f9525c23924751a3ca51").constData()),
                            _mm_loadu_si128((__m128i*) QByteArray::fromHex("324f585c6ffc1359ab371565d6c45f93").constData()),
                            _mm_loadu_si128((__m128i*) QByteArray::fromHex("ca7dd446af4aa70cc3c0cd5abba6aa1c").constData()),
                            _mm_loadu_si128((__m128i*) QByteArray::fromHex("1590df9b2eb6768289e57d56274c8570").constData())
    };
    for(i=0; i<a.length()/16; i++){
        tmp1 = _mm_loadu_si128(&((__m128i*)a.constData())[i]);
        tmp1 = _mm_shuffle_epi8(tmp1, BSWAP_MASK);
        X = _mm_xor_si128(X, tmp1);
        X = gfmul_reflected(X, H);
        //qInfo() << "X"<<(i+1)<<"  high :|"<<print128_hex_lanes(X)<<"|: low";
        assert_m128(_mm_shuffle_epi8(X, BSWAP_MASK), X_asserts[i], "X"+QString::number(i+1));       // Insight: in order to make the comparion with the values which are "normal" loaded in with QByteArray::fromHex, we need to convert the computed values back into the "normal" world domain (or outright, perform the assertion before the transformation into weird world happens)
    }
    int next_x = i;
    if(a.length() % 16){
        __m128i last_block = _mm_setzero_si128();
        for(j=0; j<a.length() % 16; j++){
            ((unsigned char*) &last_block)[j] = a[i*16 + j];
        }
        tmp1 = last_block;
        tmp1 = _mm_shuffle_epi8(tmp1, BSWAP_MASK);
        X = _mm_xor_si128(X, tmp1);
        X = gfmul_reflected(X, H);
        //qInfo() << "X"<<(i+1)<<"l high :|"<<print128_hex_lanes(X)<<"|: low";
        assert_m128(_mm_shuffle_epi8(X, BSWAP_MASK), X_asserts[i], "X"+QString::number(i+1));
        next_x = i+1;
    }

    ctr1 = j0;
    __m128i C_pre_asserts[4] = {_mm_loadu_si128((__m128i*) QByteArray::fromHex("9bb22ce7d9f372c1ee2b28722b25f206").constData()),
                            _mm_loadu_si128((__m128i*) QByteArray::fromHex("650d887c3936533a1b8d4e1ea39d2b5c").constData()),
                            _mm_loadu_si128((__m128i*) QByteArray::fromHex("3de91827c10e9a4f5240647ee5221f20").constData()),
                            _mm_loadu_si128((__m128i*) QByteArray::fromHex("aac9e6ccc0074ac0873b9ba85d908bd0").constData())
    };
    for(i=0; i<p.length() / 16; i++) {      // we traverse all **full** 16byte blocks of p
        ctr1 = _mm_shuffle_epi8(ctr1, BSWAP_EPI64_MASK);
        ctr1 = _mm_add_epi32(ctr1, ONE);
        ctr1 = _mm_shuffle_epi8(ctr1, BSWAP_EPI64_MASK);

        tmp1 = singleAESBlock(ctr1, KEY, gcm_n_rounds);     // encrypt (j0 + 1) -> tmp1
        assert_m128(tmp1, C_pre_asserts[i], "C_pre"+QString::number(i+1));

        tmp1 = _mm_xor_si128(tmp1, _mm_loadu_si128(&((__m128i*)p.constData())[i]));     // xor tmp1 with p block -> c block

        _mm_storeu_si128(&((__m128i*)out)[i], tmp1);        // store c block in out

        tmp1 = _mm_shuffle_epi8(tmp1, BSWAP_MASK);    // bring c block to reverse world

        X = _mm_xor_si128(X, tmp1);
        X = gfmul_reflected(X, H);
        assert_m128(_mm_shuffle_epi8(X, BSWAP_MASK), X_asserts[next_x++], "X"+QString::number(next_x+1));
    }
    if(p.length() % 16){
        ctr1 = _mm_shuffle_epi8(ctr1, BSWAP_EPI64_MASK);
        ctr1 = _mm_add_epi32(ctr1, ONE);
        ctr1 = _mm_shuffle_epi8(ctr1, BSWAP_EPI64_MASK);

        tmp1 = singleAESBlock(ctr1, KEY, gcm_n_rounds);     // encrypt (j0 + 1) -> tmp1
        assert_m128(tmp1, C_pre_asserts[i], "C_pre"+QString::number(i+1));

        tmp1 = _mm_xor_si128(tmp1, _mm_loadu_si128(&((__m128i*)p.constData())[i]));
        __m128i last_block = tmp1;
        for(j=0; j<p.length() % 16; j++){
            out[i*16 + j] = ((unsigned char*) &last_block)[j];
        }
        for(j; j<16; j++){
            ((unsigned char*) &last_block)[j] = 0;
        }
        tmp1 = last_block;
        tmp1 = _mm_shuffle_epi8(tmp1, BSWAP_MASK);
        X = _mm_xor_si128(X, tmp1);
        X = gfmul_reflected(X, H);
        //qInfo() << "X"<<(i+1)<<"l high :|"<<print128_hex_lanes(X)<<"|: low";
        assert_m128(_mm_shuffle_epi8(X, BSWAP_MASK), X_asserts[next_x++], "X"+QString::number(next_x+1));
    }

    tmp1 = _mm_insert_epi64(tmp1, p.length()*8, 0);
    tmp1 = _mm_insert_epi64(tmp1, a.length()*8, 1);

    X = _mm_xor_si128(X, tmp1);
    X = gfmul_reflected(X, H);
    assert_m128(_mm_shuffle_epi8(X, BSWAP_MASK), _mm_loadu_si128((__m128i*)QByteArray::fromHex("698e57f70e6ecc7fd9463b7260a9ae5f").constData()), "GHASH");
    X = _mm_shuffle_epi8(X, BSWAP_MASK);
    T = _mm_xor_si128(X, T);
    assert_m128(T, _mm_loadu_si128((__m128i*)QByteArray::fromHex("5bc94fbc3221a5db94fae95ae7121a47").constData()), "T");
    _mm_storeu_si128((__m128i*)tag, T);

    QByteArray C_assert = QByteArray::fromHex("42831ec2217774244b7221b784d0d49ce3aa212f2c02a4e035c17e2329aca12e21d514b25466931c7d8f6a5aac84aa051ba30b396a0aac973d58e091");
    if(C_assert == QByteArray::fromRawData((const char*) out, p.length())){
        qInfo() << "Assertion ( \"C\" ): OK";
    } else{
        qInfo() << "Assertion ( \"C\" ): WRONG";
    }




    /*
    ctr1 = _mm_shuffle_epi8(j0, BSWAP_EPI64_MASK);      // We need the BSWAP_EPI64 mask to bring the counter back into LE format, such that we can apply the addition and we use epi32 such that the addition stays inside 2^32
    ctr1 = _mm_add_epi32(ctr1, ONE);
    ctr2 = _mm_add_epi32(ctr1, ONE);
    ctr3 = _mm_add_epi32(ctr2, ONE);
    ctr4 = _mm_add_epi32(ctr3, ONE);


    for(i=0; i < p.length()/16/4; i++){
        tmp1 = _mm_shuffle_epi8(ctr1, BSWAP_EPI64_MASK);     // reverse the BSWAP_EPI64_MASK
        tmp2 = _mm_shuffle_epi8(ctr2, BSWAP_EPI64_MASK);
        tmp3 = _mm_shuffle_epi8(ctr3, BSWAP_EPI64_MASK);
        tmp4 = _mm_shuffle_epi8(ctr4, BSWAP_EPI64_MASK);

        ctr1 = _mm_add_epi32(ctr1, FOUR);               // bring the CTRs to the next ones (+4 on each)
        ctr2 = _mm_add_epi32(ctr2, FOUR);
        ctr3 = _mm_add_epi32(ctr3, FOUR);
        ctr4 = _mm_add_epi32(ctr4, FOUR);
        tmp1 =_mm_xor_si128(tmp1, KEY[0]);              // apply round-key 0
        tmp2 =_mm_xor_si128(tmp2, KEY[0]);
        tmp3 =_mm_xor_si128(tmp3, KEY[0]);
        tmp4 =_mm_xor_si128(tmp4, KEY[0]);
        for(j=1; j < aesKey.nr-1; j+=2){
            tmp1 = _mm_aesenc_si128(tmp1, KEY[j]);
            tmp2 = _mm_aesenc_si128(tmp2, KEY[j]);
            tmp3 = _mm_aesenc_si128(tmp3, KEY[j]);
            tmp4 = _mm_aesenc_si128(tmp4, KEY[j]);
            tmp1 = _mm_aesenc_si128(tmp1, KEY[j+1]);
            tmp2 = _mm_aesenc_si128(tmp2, KEY[j+1]);
            tmp3 = _mm_aesenc_si128(tmp3, KEY[j+1]);
            tmp4 = _mm_aesenc_si128(tmp4, KEY[j+1]);
        }
        tmp1 = _mm_aesenc_si128(tmp1, KEY[aesKey.nr-1]);
        tmp2 = _mm_aesenc_si128(tmp2, KEY[aesKey.nr-1]);
        tmp3 = _mm_aesenc_si128(tmp3, KEY[aesKey.nr-1]);
        tmp4 = _mm_aesenc_si128(tmp4, KEY[aesKey.nr-1]);
        tmp1 =_mm_aesenclast_si128(tmp1, KEY[aesKey.nr]);
        tmp2 =_mm_aesenclast_si128(tmp2, KEY[aesKey.nr]);
        tmp3 =_mm_aesenclast_si128(tmp3, KEY[aesKey.nr]);
        tmp4 =_mm_aesenclast_si128(tmp4, KEY[aesKey.nr]);
        tmp1 = _mm_xor_si128(tmp1, _mm_loadu_si128(&((__m128i*)p.constData())[i*4+0]));
        tmp2 = _mm_xor_si128(tmp2, _mm_loadu_si128(&((__m128i*)p.constData())[i*4+1]));
        tmp3 = _mm_xor_si128(tmp3, _mm_loadu_si128(&((__m128i*)p.constData())[i*4+2]));
        tmp4 = _mm_xor_si128(tmp4, _mm_loadu_si128(&((__m128i*)p.constData())[i*4+3]));
        _mm_storeu_si128(&((__m128i*)out)[i*4+0], tmp1);
        _mm_storeu_si128(&((__m128i*)out)[i*4+1], tmp2);
        _mm_storeu_si128(&((__m128i*)out)[i*4+2], tmp3);
        _mm_storeu_si128(&((__m128i*)out)[i*4+3], tmp4);
        tmp1 = _mm_shuffle_epi8(tmp1, BSWAP_MASK);
        tmp2 = _mm_shuffle_epi8(tmp2, BSWAP_MASK);
        tmp3 = _mm_shuffle_epi8(tmp3, BSWAP_MASK);
        tmp4 = _mm_shuffle_epi8(tmp4, BSWAP_MASK);
        X = _mm_xor_si128(X, tmp1);
        X= gfmul_reflected(X, H);
        X = _mm_xor_si128(X, tmp2);
        X= gfmul_reflected(X, H);
        X = _mm_xor_si128(X, tmp3);
        X= gfmul_reflected(X, H);
        X = _mm_xor_si128(X, tmp4);
        X= gfmul_reflected(X, H);
    }
    int k;
    for(k = i*4; k < p.length()/16; k++){
        tmp1 = _mm_shuffle_epi8(ctr1, BSWAP_EPI64_MASK);
        ctr1 = _mm_add_epi32(ctr1, ONE);
        tmp1 = _mm_xor_si128(tmp1, KEY[0]);
        for(j=1; j<aesKey.nr-1; j+=2){
            tmp1 = _mm_aesenc_si128(tmp1, KEY[j]);
            tmp1 = _mm_aesenc_si128(tmp1, KEY[j+1]);
        }
        tmp1 = _mm_aesenc_si128(tmp1, KEY[aesKey.nr-1]);
        tmp1 = _mm_aesenclast_si128(tmp1, KEY[aesKey.nr]);
        tmp1 = _mm_xor_si128(tmp1, _mm_loadu_si128(&((__m128i*)p.constData())[k]));
        _mm_storeu_si128(&((__m128i*)out)[k], tmp1);
        tmp1 = _mm_shuffle_epi8(tmp1, BSWAP_MASK);
        X =_mm_xor_si128(X, tmp1);
        X = gfmul_reflected(X, H);
    }

    //If one partial block remains
    if(p.length()%16){
        tmp1 = _mm_shuffle_epi8(ctr1, BSWAP_EPI64_MASK);
        tmp1 = _mm_xor_si128(tmp1, KEY[0]);
        for(j=1; j<aesKey.nr-1; j+=2){
            tmp1 = _mm_aesenc_si128(tmp1, KEY[j]);
            tmp1 = _mm_aesenc_si128(tmp1, KEY[j+1]);
        }
        tmp1 = _mm_aesenc_si128(tmp1, KEY[aesKey.nr-1]);
        tmp1 = _mm_aesenclast_si128(tmp1, KEY[aesKey.nr]);
        tmp1 = _mm_xor_si128(tmp1, _mm_loadu_si128(&((__m128i*)p.constData())[k]));
        __m128i last_block = tmp1;
        for(j=0; j < p.length()%16; j++)
            out[k*16+j]=((unsigned char*)&last_block)[j];
        for(j; j<16; j++)
            ((unsigned char*)&last_block)[j]=0;
        tmp1 = last_block;
        tmp1 = _mm_shuffle_epi8(tmp1, BSWAP_MASK);
        X =_mm_xor_si128(X, tmp1);
        X = gfmul_reflected(X, H);
    }
    tmp1 = _mm_insert_epi64(tmp1, p.length()*8, 0);
    tmp1 = _mm_insert_epi64(tmp1, a.length()*8, 1);
    X = _mm_xor_si128(X, tmp1);
    X = gfmul_reflected(X, H);
    X = _mm_shuffle_epi8(X, BSWAP_MASK);
    T = _mm_xor_si128(X, T);
    */
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
