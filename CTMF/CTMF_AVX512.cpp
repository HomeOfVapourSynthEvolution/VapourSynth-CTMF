#ifdef CTMF_X86
#include "CTMF.h"

template<typename pixel_t>
void filterRadius2_avx512(const VSFrameRef * src, VSFrameRef * dst, const int plane, const VSAPI * vsapi) noexcept {
    using vec_t = std::conditional_t<std::is_same_v<pixel_t, uint8_t>, Vec64uc, std::conditional_t<std::is_same_v<pixel_t, uint16_t>, Vec32us, Vec16f>>;

    auto sort = [](vec_t & a, vec_t & b) noexcept {
        const auto t = a;
        a = min(t, b);
        b = max(t, b);
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
        for (auto x = 0; x < width; x += vec_t().size()) {
            vec_t a[25];
            a[0] = vec_t().load(above2 + x - 2);
            a[1] = vec_t().load(above2 + x - 1);
            a[2] = vec_t().load(above2 + x);
            a[3] = vec_t().load(above2 + x + 1);
            a[4] = vec_t().load(above2 + x + 2);
            a[5] = vec_t().load(above1 + x - 2);
            a[6] = vec_t().load(above1 + x - 1);
            a[7] = vec_t().load(above1 + x);
            a[8] = vec_t().load(above1 + x + 1);
            a[9] = vec_t().load(above1 + x + 2);
            a[10] = vec_t().load(srcp + x - 2);
            a[11] = vec_t().load(srcp + x - 1);
            a[12] = vec_t().load(srcp + x);
            a[13] = vec_t().load(srcp + x + 1);
            a[14] = vec_t().load(srcp + x + 2);
            a[15] = vec_t().load(below1 + x - 2);
            a[16] = vec_t().load(below1 + x - 1);
            a[17] = vec_t().load(below1 + x);
            a[18] = vec_t().load(below1 + x + 1);
            a[19] = vec_t().load(below1 + x + 2);
            a[20] = vec_t().load(below2 + x - 2);
            a[21] = vec_t().load(below2 + x - 1);
            a[22] = vec_t().load(below2 + x);
            a[23] = vec_t().load(below2 + x + 1);
            a[24] = vec_t().load(below2 + x + 2);

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
            a[17] = max(a[8], a[17]);
            a[9] = max(a[0], a[9]);
            sort(a[10], a[19]); sort(a[1], a[19]); sort(a[1], a[10]);
            sort(a[11], a[20]); sort(a[2], a[20]); sort(a[12], a[21]);
            a[11] = max(a[2], a[11]);
            sort(a[3], a[21]); sort(a[3], a[12]); sort(a[13], a[22]);
            a[4] = min(a[4], a[22]);
            sort(a[4], a[13]); sort(a[14], a[23]);
            sort(a[5], a[23]); sort(a[5], a[14]); sort(a[15], a[24]);
            a[6] = min(a[6], a[24]);
            sort(a[6], a[15]);
            a[7] = min(a[7], a[16]);
            a[7] = min(a[7], a[19]);
            a[13] = min(a[13], a[21]);
            a[15] = min(a[15], a[23]);
            a[7] = min(a[7], a[13]);
            a[7] = min(a[7], a[15]);
            a[9] = max(a[1], a[9]);
            a[11] = max(a[3], a[11]);
            a[17] = max(a[5], a[17]);
            a[17] = max(a[11], a[17]);
            a[17] = max(a[9], a[17]);
            sort(a[4], a[10]);
            sort(a[6], a[12]); sort(a[7], a[14]); sort(a[4], a[6]);
            a[7] = max(a[4], a[7]);
            sort(a[12], a[14]);
            a[10] = min(a[10], a[14]);
            sort(a[6], a[7]); sort(a[10], a[12]); sort(a[6], a[10]);
            a[17] = max(a[6], a[17]);
            sort(a[12], a[17]);
            a[7] = min(a[7], a[17]);
            sort(a[7], a[10]); sort(a[12], a[18]);
            a[12] = max(a[7], a[12]);
            a[10] = min(a[10], a[18]);
            sort(a[12], a[20]);
            a[10] = min(a[10], a[20]);
            a[12] = max(a[10], a[12]);

            a[12].store_nt(dstp + x);
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
void ctmfHelper_avx512(const void * _srcp, void * _dstp, const CTMFData * const VS_RESTRICT d,
                       const int width, const int height, const int stride, const bool padLeft, const bool padRight) noexcept {
    /**
     * Adds histograms x and y and stores the result in y.
     */
    auto histogramAdd = [](const uint16_t * _x, uint16_t * _y) noexcept {
        for (auto i = 0; i < bins; i += Vec32us().size()) {
            const auto x = Vec32us().load_a(_x + i);
            const auto y = Vec32us().load_a(_y + i);
            (y + x).store_a(_y + i);
        }
    };

    /**
     * Subtracts histogram x from y and stores the result in y.
     */
    auto histogramSub = [](const uint16_t * _x, uint16_t * _y) noexcept {
        for (auto i = 0; i < bins; i += Vec32us().size()) {
            const auto x = Vec32us().load_a(_x + i);
            const auto y = Vec32us().load_a(_y + i);
            (y - x).store_a(_y + i);
        }
    };

    auto histogramMulAdd = [](const uint16_t a, const uint16_t * _x, uint16_t * _y) noexcept {
        for (auto i = 0; i < bins; i += Vec32us().size()) {
            const auto x = Vec32us().load_a(_x + i);
            const auto y = Vec32us().load_a(_y + i);
            (y + a * x).store_a(_y + i);
        }
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

template void filterRadius2_avx512<uint8_t>(const VSFrameRef * src, VSFrameRef * dst, const int plane, const VSAPI * vsapi) noexcept;
template void filterRadius2_avx512<uint16_t>(const VSFrameRef * src, VSFrameRef * dst, const int plane, const VSAPI * vsapi) noexcept;
template void filterRadius2_avx512<float>(const VSFrameRef * src, VSFrameRef * dst, const int plane, const VSAPI * vsapi) noexcept;

template void ctmfHelper_avx512<uint8_t, 16>(const void * srcp, void * dstp, const CTMFData * const VS_RESTRICT d, const int width, const int height, const int stride, const bool padLeft, const bool padRight) noexcept;
template void ctmfHelper_avx512<uint16_t, 32>(const void * srcp, void * dstp, const CTMFData * const VS_RESTRICT d, const int width, const int height, const int stride, const bool padLeft, const bool padRight) noexcept;
template void ctmfHelper_avx512<uint16_t, 64>(const void * srcp, void * dstp, const CTMFData * const VS_RESTRICT d, const int width, const int height, const int stride, const bool padLeft, const bool padRight) noexcept;
template void ctmfHelper_avx512<uint16_t, 128>(const void * srcp, void * dstp, const CTMFData * const VS_RESTRICT d, const int width, const int height, const int stride, const bool padLeft, const bool padRight) noexcept;
template void ctmfHelper_avx512<uint16_t, 256>(const void * srcp, void * dstp, const CTMFData * const VS_RESTRICT d, const int width, const int height, const int stride, const bool padLeft, const bool padRight) noexcept;
#endif
