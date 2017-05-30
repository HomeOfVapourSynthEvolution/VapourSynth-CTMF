#pragma once

#include <algorithm>

#include <VapourSynth.h>
#include <VSHelper.h>

#ifdef VS_TARGET_CPU_X86
#include "vectorclass/vectorclass.h"
#endif

struct CTMFData {
    VSNodeRef * node;
    const VSVideoInfo * vi;
    int radius, memsize;
    bool process[3];
    uint16_t bins, shiftRight, modulo, t;
    int stripeSize[3];
};

template<uint16_t bins>
struct alignas(32) Histogram {
    uint16_t coarse[bins];
    uint16_t fine[bins][bins];
};
