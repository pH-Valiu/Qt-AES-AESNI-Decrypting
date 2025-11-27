#ifndef GCM_H
#define GCM_H

#include "wmmintrin.h"
#include "immintrin.h"  // for _xgetbv
#include "qglobal.h"    // for Q_OS_X macros
#include <QString>
#include <QByteArray>


// ------------------------------
// AES-GCM methods
// ------------------------------

struct GCM_OUT{
    QByteArray c;
    QByteArray t;

    GCM_OUT() = default;
    GCM_OUT(QByteArray c_, QByteArray t_)
        : c(std::move(c_)), t(std::move(t_)){}
};

/**
 * @brief encrypt
 * @param key
 * @param iv
 * @param aad
 * @param p
 * @return Will return an empty GCM_OUT struct if parameters do not adhere to required lengths
 */
struct GCM_OUT encrypt(const QByteArray& key, const QByteArray& iv, const QByteArray& aad, const QByteArray& p);
QByteArray decrypt(const QByteArray& key, const QByteArray& iv, const QByteArray& aad, const QByteArray& c, const QByteArray& t);
QByteArray decrypt(const QByteArray& key, const QByteArray& iv, const QByteArray& aad, const struct GCM_OUT& gcm_out);
bool authenticate(const QByteArray& aad, const QByteArray& c, const QByteArray& t);
bool authenticate(const QByteArray& aad, const struct GCM_OUT& gcm_out);
QByteArray random_iv();

void singleAESBlock(const unsigned char* in, unsigned char* out, unsigned char length, const unsigned char* key, int number_of_rounds);

// ---------------------------------
// Future Abstraction of constants (AES-128)
// ---------------------------------
static constexpr quint8 gcm_keyLength = 16;
static constexpr quint8 gcm_keyLengthBits = 128;
static constexpr quint8 gcm_expKeyLength = 176;
static constexpr quint8 gcm_nRounds = 10;

// ----------------------------------
// AES-GCM data lengths requirements
// ----------------------------------
static constexpr quint64 MAX_PLAIN_LEN  = 549755813632ULL;          // 2^39 - 256
static constexpr quint64 MAX_AAD_LEN    = 18446744073709551615ULL;  // 2^64 - 1
static constexpr quint64 MAX_IV_LEN     = 18446744073709551615ULL;  // 2^64 - 1



// -------------------------------
// OS and CPU architecture checks
// -------------------------------

bool osSupportAVX() {
    unsigned long long xcr0 = _xgetbv(0);
    return (xcr0 & 0x6) == 0x6;      // OS AVX2 support
}
#ifdef Q_OS_WIN
#include <intrin.h>
static void cpuid(int out[], int x){
    __cpuidex(out, x, 0);
}
bool hasAES() {
    int info[4];
    cpuid(info, 1);
    return info[2] & 0x2000000 ;  // AES-NI
}
bool hasPCLMUL() {
    int info[4];
    cpuid(info, 1);
    return info[2] & 0x2;   // PCLMULQDQ
}
bool hasAVX2() {
    int info[4];
    cpuid(info, 7);
    return (info[1] & 0x20) && osSupportAVX();   // AVX2 bit
}

#elif defined(Q_OS_LINUX)
bool hasAES()       { return __builtin_cpu_supports("aes"); }
bool hasPCLMUL()    { return __builtin_cpu_supports("pclmul"); }
bool hasAVX2()      { return __builtin_cpu_supports("avx2") && osSupportAVX(); }
#else
#error "We do not support that OS version yet
#endif

#endif // GCM_H
