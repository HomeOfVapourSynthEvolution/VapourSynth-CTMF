#pragma once

#include <algorithm>
#include <memory>
#include <thread>
#include <type_traits>
#include <unordered_map>

#include <VapourSynth.h>
#include <VSHelper.h>

#ifdef CTMF_X86
#include "VCL2/vectorclass.h"
#endif

using unique_uint16 = std::unique_ptr<uint16_t[], decltype(&vs_aligned_free)>;

struct CTMFData final {
    VSNodeRef * node;
    const VSVideoInfo * vi;
    int radius;
    bool process[3];
    bool isRadius2;
    int stripeSize[3];
    uint16_t bins, shift, mask, t;
    std::unordered_map<std::thread::id, unique_uint16> hCoarse, hFine;
    void (*filter)(const VSFrameRef * src, VSFrameRef * pad, VSFrameRef * dst, const CTMFData * const VS_RESTRICT d, const VSAPI * vsapi) noexcept;
    void (*copyPad)(const VSFrameRef * src, VSFrameRef * dst, const int plane, const VSAPI * vsapi) noexcept;
    void (*filterRadius2)(const VSFrameRef * src, VSFrameRef * dst, const int plane, const VSAPI * vsapi) noexcept;
    void (*ctmfHelper)(const void * srcp, void * dstp, const CTMFData * const VS_RESTRICT d, const int width, const int height, const int stride, const bool padLeft, const bool padRight) noexcept;
};

template<uint16_t bins>
struct alignas(64) Histogram final {
    uint16_t coarse[bins];
    uint16_t fine[bins][bins];
};
