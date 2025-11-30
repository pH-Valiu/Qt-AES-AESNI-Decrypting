#ifndef GCM_H
#define GCM_H

#include "wmmintrin.h"
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
void gcm_test();

void singleAESBlock(const unsigned char* in, unsigned char* out, unsigned char length, const unsigned char* key, int number_of_rounds);
__m128i singleAESBlock(const __m128i& in, const __m128i* const key, int number_of_rounds);

// ---------------------------------
// Future Abstraction of constants (AES-128)
// ---------------------------------
static constexpr quint8 gcm_keyLength = 16;
static constexpr quint8 gcm_keyLengthBits = 128;
static constexpr quint8 gcm_expKeyLength = 176;
static constexpr quint8 gcm_n_rounds = 10;

// ----------------------------------
// AES-GCM data lengths requirements
// ----------------------------------
static constexpr quint64 MAX_PLAIN_LEN  = 549755813632ULL;          // 2^39 - 256
static constexpr quint64 MAX_AAD_LEN    = 18446744073709551615ULL;  // 2^64 - 1
static constexpr quint64 MAX_IV_LEN     = 18446744073709551615ULL;  // 2^64 - 1

// ----------------------------
// Further constant values
// ----------------------------
static __m128i BSWAP_MASK = _mm_set_epi8(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);
static __m128i BSWAP_EPI64_MASK = _mm_set_epi8(8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7);
static __m128i ONE = _mm_set_epi32(0, 1, 0, 0);
static __m128i FOUR = _mm_set_epi32(0, 4, 0, 0);


#endif // GCM_H
