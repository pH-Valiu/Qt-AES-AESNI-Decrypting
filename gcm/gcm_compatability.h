#ifndef GCM_COMPATABILITY_H
#define GCM_COMPATABILITY_H
#include "immintrin.h"  // for _xgetbv
#include "qglobal.h"    // for Q_OS_X macros


// -------------------------------
// OS and CPU architecture checks
// -------------------------------

namespace AES_GCM_Compatability{

bool osSupportAVX() {
    unsigned long long xcr0 = _xgetbv(0);
    return (xcr0 & 0x6) == 0x6;      // OS AVX2 support
}
#ifdef Q_OS_WIN
#include <intrin.h>
static void cpu_id(int out[], int x){
    __cpuidex(out, x, 0);
}
bool hasAES() {
    int info[4];
    cpu_id(info, 1);
    return info[2] & 0x2000000 ;  // AES-NI
}
bool hasPCLMUL() {
    int info[4];
    cpu_id(info, 1);
    return info[2] & 0x2;   // PCLMULQDQ
}
bool hasAVX2() {
    int info[4];
    cpu_id(info, 7);
    return (info[1] & 0x20) && osSupportAVX();   // AVX2 bit
}

#elif defined(Q_OS_LINUX)
bool hasAES()       { return __builtin_cpu_supports("aes"); }
bool hasPCLMUL()    { return __builtin_cpu_supports("pclmul"); }
bool hasAVX2()      { return __builtin_cpu_supports("avx2") && osSupportAVX(); }
#else
#error "We do not support that OS version yet"
#endif

}

#endif // GCM_COMPATABILITY_H
