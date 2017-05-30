warning_flags = -Wall -Wextra -Wshadow -Wno-unused-parameter
common_cflags = -O3 -ffast-math $(MFLAGS) $(warning_flags)
AM_CXXFLAGS = -std=c++14 $(common_cflags)

AM_CPPFLAGS = $(VapourSynth_CFLAGS)

lib_LTLIBRARIES = libctmf.la

libctmf_la_SOURCES = CTMF/CTMF.cpp \
                     CTMF/vectorclass/instrset.h \
                     CTMF/vectorclass/instrset_detect.cpp

if VS_TARGET_CPU_X86
libctmf_la_SOURCES += CTMF/CTMF_SSE2.cpp \
                      CTMF/vectorclass/vectorclass.h \
                      CTMF/vectorclass/vectorf128.h \
                      CTMF/vectorclass/vectorf256.h \
                      CTMF/vectorclass/vectorf256e.h \
                      CTMF/vectorclass/vectori128.h \
                      CTMF/vectorclass/vectori256.h \
                      CTMF/vectorclass/vectori256e.h


noinst_LTLIBRARIES = libavx2.la


libavx2_la_SOURCES = CTMF/CTMF_AVX2.cpp

libavx2_la_CXXFLAGS = $(AM_CXXFLAGS) -mfma -mavx2


libctmf_la_LIBADD = libavx2.la
endif

libctmf_la_LDFLAGS = -no-undefined -avoid-version $(PLUGINLDFLAGS)
