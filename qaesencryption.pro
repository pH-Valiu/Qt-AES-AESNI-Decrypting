QT += core testlib
QT -= gui

CONFIG += c++11

TARGET = QAESEncryption
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which as been marked deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

DEFINES += \
    USE_INTEL_AES_IF_AVAILABLE \
    AVX512_SUPPORT

QMAKE_CXXFLAGS += \
    -maes \             # enable AES-NI instruction set
    -mpclmul \          # enable Carry-Less multiply
    -mavx2 \            # enable AVX2 instruction set (includes AVX, AVX2, SSE, SSE2, SSE3, SSSE3, SSE4.1, SSE4.2)
    -mavx512f \
    -mvpclmulqdq \
    -mavx512bw \
    #-march=native \
    -msse4.1 \      # for _mm_insert_epi32
    -mssse3 \       # for _mm_shuffle_epi8
    -msse2 \        # for _mm_add_epi32 and maybe other
    -msse \         # for default compliancy for the other sse instruction sets
    #-masm=intel    ChatGPT says to not include this as it would only be needed for inline assembler but not for external files

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

HEADERS += \
    aesni/aesni-key-init.h \
    gcm/benchmarkutil.h \
    gcm/gcm.h \
    gcm/gcm_compatability.h \
    gcm/gfmul.h \
    qaesencryption.h \
    aesni/aesni-key-exp.h \
    aesni/aesni-enc-ecb.h \
    aesni/aesni-enc-cbc.h \
    unit_test/aestest.h

SOURCES += main.cpp \
    gcm/asm/gfmul_reflected_AVX_512.S \
    gcm/asm/gfmul_reflected_SSE.S \
    gcm/gcm.cpp \
    gcm/gfmul.cpp \
    qaesencryption.cpp \
    unit_test/aestest.cpp

RESOURCES += \
    res.qrc

DISTFILES += \
    GFMUL_Implementation.md \
