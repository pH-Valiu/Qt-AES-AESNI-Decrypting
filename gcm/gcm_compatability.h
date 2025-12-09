#ifndef GCM_COMPATABILITY_H
#define GCM_COMPATABILITY_H
#include "immintrin.h"  // for _xgetbv
#include "qglobal.h"    // for Q_OS_X macros


// -------------------------------
// OS and CPU architecture checks
// -------------------------------

namespace AES_GCM_Compatability{

inline bool osSupportAVX() {
    unsigned long long xcr0 = _xgetbv(0);
    return (xcr0 & 0x6) == 0x6;      // OS AVX2 support
}
inline bool osSupportsAVX512() {
    unsigned long long xcr0 = _xgetbv(0);
    // XMM + YMM + ZMM-hi256 + ZMM-hi16 registers
    // Bits 1,2,5,6,7 must be set.
    return (xcr0 & 0xe6) == 0xe6;
}


#ifdef Q_OS_WIN
#include <intrin.h>
inline void cpu_id(int out[], int x){
    __cpuidex(out, x, 0);
}

// --- SSE family ---
inline bool hasSSE2()  { int a[4]; cpu_id(a, 1); return (a[3] & 0x4000000); }
inline bool hasSSE3()  { int a[4]; cpu_id(a, 1); return (a[2] & 0x1);  }
inline bool hasSSSE3() { int a[4]; cpu_id(a, 1); return (a[2] & 0x2);  }
inline bool hasSSE41() { int a[4]; cpu_id(a, 1); return (a[2] & 0x80000); }
inline bool hasSSE42() { int a[4]; cpu_id(a, 1); return (a[2] & 0x100000); }

inline bool hasAES() {int info[4]; cpu_id(info, 1); return info[2] & 0x2000000; }  // AES-NI
inline bool hasPCLMUL() {int info[4]; cpu_id(info, 1); return info[2] & 0x2; }   // PCLMULQDQ
inline bool hasAVX2() {int info[4]; cpu_id(info, 7); return (info[1] & 0x20) && osSupportAVX(); }   // AVX2 bit

// --- AVX-512: requires AVX-512F plus OS support ---
inline bool hasAVX512F()   {int a[4]; cpu_id(a, 7); return (a[1] & (1 << 16)) && osSupportsAVX512(); } //AVX-512-F
inline bool hasAVX512VL()  {int a[4]; cpu_id(a, 7); return (a[1] & (1 << 31)) && hasAVX512F(); } // AVX-512-VL
inline bool hasAVX512DQ()  {int a[4]; cpu_id(a, 7); return (a[1] & (1 << 17)) && hasAVX512F(); } // AVX-512-DQ
inline bool hasAVX512BW()  {int a[4]; cpu_id(a, 7); return (a[1] & (1 << 30)) && hasAVX512F(); } // AVX-512-BW

inline bool hasAVX512()    {return hasAVX512F() && hasAVX512VL() && hasAVX512BW() && hasAVX512DQ(); }

#elif defined(Q_OS_LINUX)
inline bool hasSSE2()      { return __builtin_cpu_supports("sse2"); }
inline bool hasSSE3()      { return __builtin_cpu_supports("sse3"); }
inline bool hasSSSE3()     { return __builtin_cpu_supports("ssse3"); }
inline bool hasSSE41()     { return __builtin_cpu_supports("sse4.1"); }
inline bool hasSSE42()     { return __builtin_cpu_supports("sse4.2"); }
inline bool hasAES()       { return __builtin_cpu_supports("aes"); }
inline bool hasPCLMUL()    { return __builtin_cpu_supports("pclmul"); }
inline bool hasAVX2()      { return __builtin_cpu_supports("avx2") && osSupportAVX(); }
inline bool hasAVX512F()   { return __builtin_cpu_supports("avx512f")  && osSupportsAVX512(); }
inline bool hasAVX512VL()  { return __builtin_cpu_supports("avx512vl") && hasAVX512F(); }
inline bool hasAVX512DQ()  { return __builtin_cpu_supports("avx512dq") && hasAVX512F(); }
inline bool hasAVX512BW()  { return __builtin_cpu_supports("avx512bw") && hasAVX512F(); }
inline bool hasAVX512()    { return hasAVX512F() && hasAVX512VL() && hasAVX512BW() && hasAVX512DQ(); }

#else
#error "We do not support that OS version yet"
#endif

}

#endif // GCM_COMPATABILITY_H
