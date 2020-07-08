/*
 * VapourSynth port by HolyWu
 *
 * ctmf.c - Constant-time median filtering
 * Copyright (C) 2006  Simon Perreault
 *
 * Reference: S. Perreault and P. Hébert, "Median Filtering in Constant Time",
 * IEEE Transactions on Image Processing, September 2007.
 *
 * This program has been obtained from http://nomis80.org/ctmf.html. No patent
 * covers this program, although it is subject to the following license:
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 * Contact:
 *  Laboratoire de vision et systèmes numériques
 *  Pavillon Adrien-Pouliot
 *  Université Laval
 *  Sainte-Foy, Québec, Canada
 *  G1K 7P4
 *
 *  perreaul@gel.ulaval.ca
 */

#include <cmath>
#include <string>

#include "CTMF.h"

using namespace std::literals;

#ifdef CTMF_X86
template<typename pixel_t> extern void filterRadius2_sse2(const VSFrameRef * src, VSFrameRef * dst, const int plane, const VSAPI * vsapi) noexcept;
template<typename pixel_t> extern void filterRadius2_avx2(const VSFrameRef * src, VSFrameRef * dst, const int plane, const VSAPI * vsapi) noexcept;
template<typename pixel_t> extern void filterRadius2_avx512(const VSFrameRef * src, VSFrameRef * dst, const int plane, const VSAPI * vsapi) noexcept;

template<typename pixel_t, uint16_t bins> extern void ctmfHelper_sse2(const void * srcp, void * dstp, const CTMFData * const VS_RESTRICT d, const int width, const int height, const int stride, const bool padLeft, const bool padRight) noexcept;
template<typename pixel_t, uint16_t bins> extern void ctmfHelper_avx2(const void * srcp, void * dstp, const CTMFData * const VS_RESTRICT d, const int width, const int height, const int stride, const bool padLeft, const bool padRight) noexcept;
template<typename pixel_t, uint16_t bins> extern void ctmfHelper_avx512(const void * srcp, void * dstp, const CTMFData * const VS_RESTRICT d, const int width, const int height, const int stride, const bool padLeft, const bool padRight) noexcept;
#endif

template<typename pixel_t>
static auto copyPad(const VSFrameRef * src, VSFrameRef * dst, const int plane, const VSAPI * vsapi) noexcept {
    const auto width = vsapi->getFrameWidth(src, plane);
    const auto height = vsapi->getFrameHeight(src, plane);
    const int stride = vsapi->getStride(dst, 0) / sizeof(pixel_t);
    auto srcp = reinterpret_cast<const pixel_t *>(vsapi->getReadPtr(src, plane));
    auto dstp = reinterpret_cast<pixel_t *>(vsapi->getWritePtr(dst, 0)) + stride * 2 + 2;

    vs_bitblt(dstp, vsapi->getStride(dst, 0), srcp, vsapi->getStride(src, plane), width * sizeof(pixel_t), height);

    for (auto y = 0; y < height; y++) {
        dstp[-1] = dstp[-2] = 0;
        dstp[width] = dstp[width + 1] = 0;

        dstp += stride;
    }

    dstp = reinterpret_cast<pixel_t *>(vsapi->getWritePtr(dst, 0));
    memset(dstp, 0, vsapi->getStride(dst, 0) * 2);
    dstp += stride * (vsapi->getFrameHeight(dst, 0) - 2);
    memset(dstp, 0, vsapi->getStride(dst, 0) * 2);
}

template<typename pixel_t>
static void filterRadius2_c(const VSFrameRef * src, VSFrameRef * dst, const int plane, const VSAPI * vsapi) noexcept {
    auto sort = [](pixel_t & a, pixel_t & b) noexcept {
        const auto t = a;
        a = std::min(t, b);
        b = std::max(t, b);
    };

    const auto width = vsapi->getFrameWidth(dst, plane);
    const auto height = vsapi->getFrameHeight(dst, plane);
    const int srcStride = vsapi->getStride(src, 0) / sizeof(pixel_t);
    const int dstStride = vsapi->getStride(dst, plane) / sizeof(pixel_t);
    auto srcp = reinterpret_cast<const pixel_t *>(vsapi->getReadPtr(src, 0)) + srcStride * 2 + 2;
    auto dstp = reinterpret_cast<pixel_t *>(vsapi->getWritePtr(dst, plane));

    auto above2 = srcp - srcStride * 2;
    auto above1 = srcp - srcStride;
    auto below1 = srcp + srcStride;
    auto below2 = srcp + srcStride * 2;

    for (auto y = 0; y < height; y++) {
        for (auto x = 0; x < width; x++) {
            pixel_t a[25];
            a[0] = above2[x - 2]; a[1] = above2[x - 1]; a[2] = above2[x]; a[3] = above2[x + 1]; a[4] = above2[x + 2];
            a[5] = above1[x - 2]; a[6] = above1[x - 1]; a[7] = above1[x]; a[8] = above1[x + 1]; a[9] = above1[x + 2];
            a[10] = srcp[x - 2]; a[11] = srcp[x - 1]; a[12] = srcp[x]; a[13] = srcp[x + 1]; a[14] = srcp[x + 2];
            a[15] = below1[x - 2]; a[16] = below1[x - 1]; a[17] = below1[x]; a[18] = below1[x + 1]; a[19] = below1[x + 2];
            a[20] = below2[x - 2]; a[21] = below2[x - 1]; a[22] = below2[x]; a[23] = below2[x + 1]; a[24] = below2[x + 2];

            sort(a[0], a[1]); sort(a[3], a[4]); sort(a[2], a[4]);
            sort(a[2], a[3]); sort(a[6], a[7]); sort(a[5], a[7]);
            sort(a[5], a[6]); sort(a[9], a[10]); sort(a[8], a[10]);
            sort(a[8], a[9]); sort(a[12], a[13]); sort(a[11], a[13]);
            sort(a[11], a[12]); sort(a[15], a[16]); sort(a[14], a[16]);
            sort(a[14], a[15]); sort(a[18], a[19]); sort(a[17], a[19]);
            sort(a[17], a[18]); sort(a[21], a[22]); sort(a[20], a[22]);
            sort(a[20], a[21]); sort(a[23], a[24]); sort(a[2], a[5]);
            sort(a[3], a[6]); sort(a[0], a[6]); sort(a[0], a[3]);
            sort(a[4], a[7]); sort(a[1], a[7]); sort(a[1], a[4]);
            sort(a[11], a[14]); sort(a[8], a[14]); sort(a[8], a[11]);
            sort(a[12], a[15]); sort(a[9], a[15]); sort(a[9], a[12]);
            sort(a[13], a[16]); sort(a[10], a[16]); sort(a[10], a[13]);
            sort(a[20], a[23]); sort(a[17], a[23]); sort(a[17], a[20]);
            sort(a[21], a[24]); sort(a[18], a[24]); sort(a[18], a[21]);
            sort(a[19], a[22]); sort(a[9], a[18]); sort(a[0], a[18]);
            a[17] = std::max(a[8], a[17]);
            a[9] = std::max(a[0], a[9]);
            sort(a[10], a[19]); sort(a[1], a[19]); sort(a[1], a[10]);
            sort(a[11], a[20]); sort(a[2], a[20]); sort(a[12], a[21]);
            a[11] = std::max(a[2], a[11]);
            sort(a[3], a[21]); sort(a[3], a[12]); sort(a[13], a[22]);
            a[4] = std::min(a[4], a[22]);
            sort(a[4], a[13]); sort(a[14], a[23]);
            sort(a[5], a[23]); sort(a[5], a[14]); sort(a[15], a[24]);
            a[6] = std::min(a[6], a[24]);
            sort(a[6], a[15]);
            a[7] = std::min(a[7], a[16]);
            a[7] = std::min(a[7], a[19]);
            a[13] = std::min(a[13], a[21]);
            a[15] = std::min(a[15], a[23]);
            a[7] = std::min(a[7], a[13]);
            a[7] = std::min(a[7], a[15]);
            a[9] = std::max(a[1], a[9]);
            a[11] = std::max(a[3], a[11]);
            a[17] = std::max(a[5], a[17]);
            a[17] = std::max(a[11], a[17]);
            a[17] = std::max(a[9], a[17]);
            sort(a[4], a[10]);
            sort(a[6], a[12]); sort(a[7], a[14]); sort(a[4], a[6]);
            a[7] = std::max(a[4], a[7]);
            sort(a[12], a[14]);
            a[10] = std::min(a[10], a[14]);
            sort(a[6], a[7]); sort(a[10], a[12]); sort(a[6], a[10]);
            a[17] = std::max(a[6], a[17]);
            sort(a[12], a[17]);
            a[7] = std::min(a[7], a[17]);
            sort(a[7], a[10]); sort(a[12], a[18]);
            a[12] = std::max(a[7], a[12]);
            a[10] = std::min(a[10], a[18]);
            sort(a[12], a[20]);
            a[10] = std::min(a[10], a[20]);
            a[12] = std::max(a[10], a[12]);

            dstp[x] = a[12];
        }

        above2 += srcStride;
        above1 += srcStride;
        srcp += srcStride;
        below1 += srcStride;
        below2 += srcStride;
        dstp += dstStride;
    }
}

template<typename pixel_t, uint16_t bins>
static void ctmfHelper_c(const void * _srcp, void * _dstp, const CTMFData * const VS_RESTRICT d,
                         const int width, const int height, const int stride, const bool padLeft, const bool padRight) noexcept {
    /**
     * Adds histograms x and y and stores the result in y.
     */
    auto histogramAdd = [](const uint16_t * x, uint16_t * VS_RESTRICT y) noexcept {
        for (auto i = 0; i < bins; i++)
            y[i] += x[i];
    };

    /**
     * Subtracts histogram x from y and stores the result in y.
     */
    auto histogramSub = [](const uint16_t * x, uint16_t * VS_RESTRICT y) noexcept {
        for (auto i = 0; i < bins; i++)
            y[i] -= x[i];
    };

    auto histogramMulAdd = [](const uint16_t a, const uint16_t * x, uint16_t * VS_RESTRICT y) noexcept {
        for (auto i = 0; i < bins; i++)
            y[i] += a * x[i];
    };

    auto srcp = static_cast<const pixel_t *>(_srcp);
    auto dstp = static_cast<pixel_t *>(_dstp);

    const auto threadId = std::this_thread::get_id();
    auto hCoarse = d->hCoarse.at(threadId).get();
    auto hFine = d->hFine.at(threadId).get();

    const pixel_t * p, * q;

    Histogram<bins> H;
    uint16_t luc[bins];

    memset(hCoarse, 0, bins * width * sizeof(uint16_t));
    memset(hFine, 0, bins * bins * width * sizeof(uint16_t));

    // First row initialization
    for (auto j = 0; j < width; j++) {
        hCoarse[bins * j + (srcp[j] >> d->shift)] += d->radius + 1;
        hFine[bins * (width * (srcp[j] >> d->shift) + j) + (srcp[j] & d->mask)] += d->radius + 1;
    }
    for (auto i = 0; i < d->radius; i++) {
        for (auto j = 0; j < width; j++) {
            hCoarse[bins * j + (srcp[stride * i + j] >> d->shift)]++;
            hFine[bins * (width * (srcp[stride * i + j] >> d->shift) + j) + (srcp[stride * i + j] & d->mask)]++;
        }
    }

    for (auto i = 0; i < height; i++) {
        // Update column histograms for entire row
        p = srcp + stride * std::max(0, i - d->radius - 1);
        q = p + width;
        for (auto j = 0; p != q; j++) {
            hCoarse[bins * j + (*p >> d->shift)]--;
            hFine[bins * (width * (*p >> d->shift) + j) + (*p & d->mask)]--;
            p++;
        }

        p = srcp + stride * std::min(height - 1, i + d->radius);
        q = p + width;
        for (auto j = 0; p != q; j++) {
            hCoarse[bins * j + (*p >> d->shift)]++;
            hFine[bins * (width * (*p >> d->shift) + j) + (*p & d->mask)]++;
            p++;
        }

        // First column initialization
        memset(&H, 0, sizeof(H));
        memset(luc, 0, sizeof(luc));
        if (padLeft)
            histogramMulAdd(d->radius, &hCoarse[0], H.coarse);
        for (auto j = 0; j < (padLeft ? d->radius : 2 * d->radius); j++)
            histogramAdd(&hCoarse[bins * j], H.coarse);
        for (auto k = 0; k < bins; k++)
            histogramMulAdd(2 * d->radius + 1, &hFine[bins * width * k], &H.fine[k][0]);

        for (auto j = padLeft ? 0 : d->radius; j < (padRight ? width : width - d->radius); j++) {
            uint16_t sum = 0, * segment;
            int k;

            histogramAdd(&hCoarse[bins * std::min(j + d->radius, width - 1)], H.coarse);

            // Find median at coarse level
            for (k = 0; k < bins; k++) {
                sum += H.coarse[k];
                if (sum > d->t) {
                    sum -= H.coarse[k];
                    break;
                }
            }

            // Update corresponding histogram segment
            if (luc[k] <= j - d->radius) {
                memset(&H.fine[k], 0, bins * sizeof(uint16_t));
                for (luc[k] = j - d->radius; luc[k] < std::min(j + d->radius + 1, width); luc[k]++)
                    histogramAdd(&hFine[bins * (width * k + luc[k])], H.fine[k]);
                if (luc[k] < j + d->radius + 1) {
                    histogramMulAdd(j + d->radius + 1 - width, &hFine[bins * (width * k + width - 1)], &H.fine[k][0]);
                    luc[k] = j + d->radius + 1;
                }
            } else {
                for (; luc[k] < j + d->radius + 1; luc[k]++) {
                    histogramSub(&hFine[bins * (width * k + std::max(luc[k] - 2 * d->radius - 1, 0))], H.fine[k]);
                    histogramAdd(&hFine[bins * (width * k + std::min(static_cast<int>(luc[k]), width - 1))], H.fine[k]);
                }
            }

            histogramSub(&hCoarse[bins * std::max(j - d->radius, 0)], H.coarse);

            // Find median in segment
            segment = H.fine[k];
            for (auto b = 0; b < bins; b++) {
                sum += segment[b];
                if (sum > d->t) {
                    dstp[stride * i + j] = bins * k + b;
                    break;
                }
            }
        }
    }
}

template<typename pixel_t>
static void filter(const VSFrameRef * src, VSFrameRef * pad, VSFrameRef * dst, const CTMFData * const VS_RESTRICT d, const VSAPI * vsapi) noexcept {
    for (auto plane = 0; plane < d->vi->format->numPlanes; plane++) {
        if (d->process[plane]) {
            if (d->isRadius2) {
                d->copyPad(src, pad, plane, vsapi);
                d->filterRadius2(pad, dst, plane, vsapi);
            } else {
                const auto width = vsapi->getFrameWidth(src, plane);
                const auto height = vsapi->getFrameHeight(src, plane);
                const int stride = vsapi->getStride(src, plane) / sizeof(pixel_t);
                auto srcp = reinterpret_cast<const pixel_t *>(vsapi->getReadPtr(src, plane));
                auto dstp = reinterpret_cast<pixel_t *>(vsapi->getWritePtr(dst, plane));

                for (auto i = 0; i < width; i += d->stripeSize[plane] - 2 * d->radius) {
                    auto stripe = d->stripeSize[plane];
                    // Make sure that the filter kernel fits into one stripe
                    if (i + d->stripeSize[plane] - 2 * d->radius >= width || width - (i + d->stripeSize[plane] - 2 * d->radius) < 2 * d->radius + 1)
                        stripe = width - i;

                    d->ctmfHelper(srcp + i, dstp + i, d, stripe, height, stride, i == 0, stripe == width - i);

                    if (stripe == width - i)
                        break;
                }
            }
        }
    }
}

static void VS_CC ctmfInit(VSMap * in, VSMap * out, void ** instanceData, VSNode * node, VSCore * core, const VSAPI * vsapi) {
    auto d = static_cast<CTMFData *>(*instanceData);
    vsapi->setVideoInfo(d->vi, 1, node);
}

static const VSFrameRef * VS_CC ctmfGetFrame(int n, int activationReason, void ** instanceData, void ** frameData, VSFrameContext * frameCtx, VSCore * core, const VSAPI * vsapi) {
    auto d = static_cast<CTMFData *>(*instanceData);

    if (activationReason == arInitial) {
        vsapi->requestFrameFilter(n, d->node, frameCtx);
    } else if (activationReason == arAllFramesReady) {
        try {
            auto threadId = std::this_thread::get_id();

            if (!d->isRadius2 && !d->hCoarse.count(threadId)) {
                auto hCoarse = vs_aligned_malloc<uint16_t>(d->bins * d->vi->width * sizeof(uint16_t), 64);
                if (!hCoarse)
                    throw "malloc failure (hCoarse)";
                d->hCoarse.emplace(threadId, unique_uint16{ hCoarse, vs_aligned_free });

                auto hFine = vs_aligned_malloc<uint16_t>(d->bins * d->bins * d->vi->width * sizeof(uint16_t), 64);
                if (!hFine)
                    throw "malloc failure (hFine)";
                d->hFine.emplace(threadId, unique_uint16{ hFine, vs_aligned_free });
            }
        } catch (const char * error) {
            vsapi->setFilterError(("CTMF: "s + error).c_str(), frameCtx);
            return nullptr;
        }

        auto src = vsapi->getFrameFilter(n, d->node, frameCtx);
        const VSFrameRef * fr[] = { d->process[0] ? nullptr : src, d->process[1] ? nullptr : src, d->process[2] ? nullptr : src };
        const int pl[] = { 0, 1, 2 };
        auto dst = vsapi->newVideoFrame2(d->vi->format, d->vi->width, d->vi->height, fr, pl, src, core);

        auto pad = [&]() noexcept -> VSFrameRef * {
            if (d->isRadius2)
                return vsapi->newVideoFrame(vsapi->registerFormat(cmGray, d->vi->format->sampleType, d->vi->format->bitsPerSample, 0, 0, core),
                                            d->vi->width + 4, d->vi->height + 4, nullptr, core);
            else
                return nullptr;
        }();

        d->filter(src, pad, dst, d, vsapi);

        vsapi->freeFrame(src);
        vsapi->freeFrame(pad);
        return dst;
    }

    return nullptr;
}

static void VS_CC ctmfFree(void * instanceData, VSCore * core, const VSAPI * vsapi) {
    auto d = static_cast<CTMFData *>(instanceData);
    vsapi->freeNode(d->node);
    delete d;
}

static void VS_CC ctmfCreate(const VSMap * in, VSMap * out, void * userData, VSCore * core, const VSAPI * vsapi) {
    auto d = std::make_unique<CTMFData>();
    int err;

    try {
        d->node = vsapi->propGetNode(in, "clip", 0, nullptr);
        d->vi = vsapi->getVideoInfo(d->node);

        if (!isConstantFormat(d->vi) ||
            (d->vi->format->sampleType == stInteger && (d->vi->format->bitsPerSample > 16 || d->vi->format->bitsPerSample & 1)) ||
            (d->vi->format->sampleType == stFloat && d->vi->format->bitsPerSample != 32))
            throw "only constant format 8, 10, 12, 14, 16 bit integer and 32 bit float input supported";

        d->radius = int64ToIntS(vsapi->propGetInt(in, "radius", 0, &err));
        if (err)
            d->radius = 2;

        auto memsize = int64ToIntS(vsapi->propGetInt(in, "memsize", 0, &err));
        if (err)
            memsize = 1048576;

        {
            const auto m = vsapi->propNumElements(in, "planes");

            for (int i = 0; i < 3; i++)
                d->process[i] = (m <= 0);

            for (int i = 0; i < m; i++) {
                const auto n = int64ToIntS(vsapi->propGetInt(in, "planes", i, nullptr));

                if (n < 0 || n >= d->vi->format->numPlanes)
                    throw "plane index out of range";

                if (d->process[n])
                    throw "plane specified twice";

                d->process[n] = true;
            }
        }

        const auto opt = int64ToIntS(vsapi->propGetInt(in, "opt", 0, &err));

        if (d->radius < 1 || d->radius > 127)
            throw "radius must be between 1 and 127 (inclusive)";

        if (d->vi->format->sampleType == stFloat && d->radius != 2)
            throw "32 bit float input is only supported for radius=2";

        if (memsize < 1024)
            throw "memsize must be greater than or equal to 1024";

        if (opt < 0 || opt > 4)
            throw "opt must be 0, 1, 2, 3, or 4";

        d->isRadius2 = (d->radius == 2);

        {
            if (d->vi->format->bytesPerSample == 1) {
                d->filter = filter<uint8_t>;
                if (d->isRadius2) {
                    d->copyPad = copyPad<uint8_t>;
                    d->filterRadius2 = filterRadius2_c<uint8_t>;
                }
            } else if (d->vi->format->bytesPerSample == 2) {
                d->filter = filter<uint16_t>;
                if (d->isRadius2) {
                    d->copyPad = copyPad<uint16_t>;
                    d->filterRadius2 = filterRadius2_c<uint16_t>;
                }
            } else {
                d->filter = filter<float>;
                if (d->isRadius2) {
                    d->copyPad = copyPad<float>;
                    d->filterRadius2 = filterRadius2_c<float>;
                }
            }

            if (!d->isRadius2) {
                if (d->vi->format->bitsPerSample == 8)
                    d->ctmfHelper = ctmfHelper_c<uint8_t, 16>;
                else if (d->vi->format->bitsPerSample == 10)
                    d->ctmfHelper = ctmfHelper_c<uint16_t, 32>;
                else if (d->vi->format->bitsPerSample == 12)
                    d->ctmfHelper = ctmfHelper_c<uint16_t, 64>;
                else if (d->vi->format->bitsPerSample == 14)
                    d->ctmfHelper = ctmfHelper_c<uint16_t, 128>;
                else if (d->vi->format->bitsPerSample == 16)
                    d->ctmfHelper = ctmfHelper_c<uint16_t, 256>;
            }

#ifdef CTMF_X86
            const int iset = instrset_detect();
            if ((opt == 0 && iset >= 10) || opt == 4) {
                if (d->isRadius2) {
                    if (d->vi->format->bytesPerSample == 1)
                        d->filterRadius2 = filterRadius2_avx512<uint8_t>;
                    else if (d->vi->format->bytesPerSample == 2)
                        d->filterRadius2 = filterRadius2_avx512<uint16_t>;
                    else
                        d->filterRadius2 = filterRadius2_avx512<float>;
                } else {
                    if (d->vi->format->bitsPerSample == 8)
                        d->ctmfHelper = ctmfHelper_avx512<uint8_t, 16>;
                    else if (d->vi->format->bitsPerSample == 10)
                        d->ctmfHelper = ctmfHelper_avx512<uint16_t, 32>;
                    else if (d->vi->format->bitsPerSample == 12)
                        d->ctmfHelper = ctmfHelper_avx512<uint16_t, 64>;
                    else if (d->vi->format->bitsPerSample == 14)
                        d->ctmfHelper = ctmfHelper_avx512<uint16_t, 128>;
                    else if (d->vi->format->bitsPerSample == 16)
                        d->ctmfHelper = ctmfHelper_avx512<uint16_t, 256>;
                }
            } else if ((opt == 0 && iset >= 8) || opt == 3) {
                if (d->isRadius2) {
                    if (d->vi->format->bytesPerSample == 1)
                        d->filterRadius2 = filterRadius2_avx2<uint8_t>;
                    else if (d->vi->format->bytesPerSample == 2)
                        d->filterRadius2 = filterRadius2_avx2<uint16_t>;
                    else
                        d->filterRadius2 = filterRadius2_avx2<float>;
                } else {
                    if (d->vi->format->bitsPerSample == 8)
                        d->ctmfHelper = ctmfHelper_avx2<uint8_t, 16>;
                    else if (d->vi->format->bitsPerSample == 10)
                        d->ctmfHelper = ctmfHelper_avx2<uint16_t, 32>;
                    else if (d->vi->format->bitsPerSample == 12)
                        d->ctmfHelper = ctmfHelper_avx2<uint16_t, 64>;
                    else if (d->vi->format->bitsPerSample == 14)
                        d->ctmfHelper = ctmfHelper_avx2<uint16_t, 128>;
                    else if (d->vi->format->bitsPerSample == 16)
                        d->ctmfHelper = ctmfHelper_avx2<uint16_t, 256>;
                }
            } else if ((opt == 0 && iset >= 2) || opt == 2) {
                if (d->isRadius2) {
                    if (d->vi->format->bytesPerSample == 1)
                        d->filterRadius2 = filterRadius2_sse2<uint8_t>;
                    else if (d->vi->format->bytesPerSample == 2)
                        d->filterRadius2 = filterRadius2_sse2<uint16_t>;
                    else
                        d->filterRadius2 = filterRadius2_sse2<float>;
                } else {
                    if (d->vi->format->bitsPerSample == 8)
                        d->ctmfHelper = ctmfHelper_sse2<uint8_t, 16>;
                    else if (d->vi->format->bitsPerSample == 10)
                        d->ctmfHelper = ctmfHelper_sse2<uint16_t, 32>;
                    else if (d->vi->format->bitsPerSample == 12)
                        d->ctmfHelper = ctmfHelper_sse2<uint16_t, 64>;
                    else if (d->vi->format->bitsPerSample == 14)
                        d->ctmfHelper = ctmfHelper_sse2<uint16_t, 128>;
                    else if (d->vi->format->bitsPerSample == 16)
                        d->ctmfHelper = ctmfHelper_sse2<uint16_t, 256>;
                }
            }
#endif
        }

        const auto histogramSize = [&]() noexcept {
            if (d->vi->format->bitsPerSample == 8)
                return sizeof(Histogram<16>);
            else if (d->vi->format->bitsPerSample == 10)
                return sizeof(Histogram<32>);
            else if (d->vi->format->bitsPerSample == 12)
                return sizeof(Histogram<64>);
            else if (d->vi->format->bitsPerSample == 14)
                return sizeof(Histogram<128>);
            else if (d->vi->format->bitsPerSample == 16)
                return sizeof(Histogram<256>);
            else
                return static_cast<size_t>(1);
        }();

        for (auto plane = 0; plane < d->vi->format->numPlanes; plane++) {
            if (d->process[plane]) {
                const auto width = d->vi->width >> (plane ? d->vi->format->subSamplingW : 0);
                const auto height = d->vi->height >> (plane ? d->vi->format->subSamplingH : 0);

                if (width < 2 * d->radius + 1)
                    throw "plane's width must be greater than or equal to 2*radius+1";

                if (height < 2 * d->radius + 1)
                    throw "plane's height must be greater than or equal to 2*radius+1";

                if (!d->isRadius2) {
                    const auto stripes = static_cast<int>(std::ceil(static_cast<float>(width - 2 * d->radius) / (memsize / histogramSize - 2 * d->radius)));
                    d->stripeSize[plane] = static_cast<int>(std::ceil(static_cast<float>(width + stripes * 2 * d->radius - 2 * d->radius) / stripes));
                }
            }
        }

        if (!d->isRadius2) {
            d->bins = 1 << (d->vi->format->bitsPerSample / 2);
            d->shift = d->vi->format->bitsPerSample / 2;
            d->mask = d->bins - 1;

            d->t = 2 * d->radius * d->radius + 2 * d->radius;

            const auto numThreads = vsapi->getCoreInfo(core)->numThreads;
            d->hCoarse.reserve(numThreads);
            d->hFine.reserve(numThreads);
        }
    } catch (const char * error) {
        vsapi->setError(out, ("CTMF: "s + error).c_str());
        vsapi->freeNode(d->node);
        return;
    }

    vsapi->createFilter(in, out, "CTMF", ctmfInit, ctmfGetFrame, ctmfFree, fmParallel, 0, d.release(), core);
}

//////////////////////////////////////////
// Init

VS_EXTERNAL_API(void) VapourSynthPluginInit(VSConfigPlugin configFunc, VSRegisterFunction registerFunc, VSPlugin * plugin) {
    configFunc("com.holywu.ctmf", "ctmf", "Constant Time Median Filtering", VAPOURSYNTH_API_VERSION, 1, plugin);
    registerFunc("CTMF",
                 "clip:clip;"
                 "radius:int:opt;"
                 "memsize:int:opt;"
                 "planes:int[]:opt;"
                 "opt:int:opt;",
                 ctmfCreate, nullptr, plugin);
}
