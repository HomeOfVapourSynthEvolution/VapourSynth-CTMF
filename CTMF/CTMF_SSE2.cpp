#ifdef VS_TARGET_CPU_X86
#include "CTMF.hpp"

static inline void histogramAdd_sse2(const uint16_t * _x, uint16_t * _y, const uint16_t bins) noexcept {
    for (uint16_t i = 0; i < bins; i += 8) {
        const Vec8us x = Vec8us().load_a(_x + i);
        const Vec8us y = Vec8us().load_a(_y + i);
        (y + x).store_a(_y + i);
    }
}

static inline void histogramSub_sse2(const uint16_t * _x, uint16_t * _y, const uint16_t bins) noexcept {
    for (uint16_t i = 0; i < bins; i += 8) {
        const Vec8us x = Vec8us().load_a(_x + i);
        const Vec8us y = Vec8us().load_a(_y + i);
        (y - x).store_a(_y + i);
    }
}

static inline void histogramMulAdd_sse2(const uint16_t a, const uint16_t * _x, uint16_t * _y, const uint16_t bins) noexcept {
    for (uint16_t i = 0; i < bins; i += 8) {
        const Vec8us x = Vec8us().load_a(_x + i);
        const Vec8us y = Vec8us().load_a(_y + i);
        (y + a * x).store_a(_y + i);
    }
}

template<typename T, uint16_t bins>
void process_sse2(const T * srcp, T * VS_RESTRICT dstp, uint16_t * VS_RESTRICT hCoarse, uint16_t * VS_RESTRICT hFine, const CTMFData * d,
                  const int width, const int height, const int stride, const bool padLeft, const bool padRight) noexcept {
    const T * p, * q;

    Histogram<bins> H;
    uint16_t luc[bins];

    memset(hCoarse, 0, bins * width * sizeof(uint16_t));
    memset(hFine, 0, bins * bins * width * sizeof(uint16_t));

    // First row initialization
    for (int j = 0; j < width; j++) {
        hCoarse[bins * j + (srcp[j] >> d->shiftRight)] += d->radius + 1;
        hFine[bins * (width * (srcp[j] >> d->shiftRight) + j) + (srcp[j] & d->modulo)] += d->radius + 1;
    }
    for (int i = 0; i < d->radius; i++) {
        for (int j = 0; j < width; j++) {
            hCoarse[bins * j + (srcp[stride * i + j] >> d->shiftRight)]++;
            hFine[bins * (width * (srcp[stride * i + j] >> d->shiftRight) + j) + (srcp[stride * i + j] & d->modulo)]++;
        }
    }

    for (int i = 0; i < height; i++) {
        // Update column histograms for entire row
        p = srcp + stride * std::max(0, i - d->radius - 1);
        q = p + width;
        for (int j = 0; p != q; j++) {
            hCoarse[bins * j + (*p >> d->shiftRight)]--;
            hFine[bins * (width * (*p >> d->shiftRight) + j) + (*p & d->modulo)]--;
            p++;
        }

        p = srcp + stride * std::min(height - 1, i + d->radius);
        q = p + width;
        for (int j = 0; p != q; j++) {
            hCoarse[bins * j + (*p >> d->shiftRight)]++;
            hFine[bins * (width * (*p >> d->shiftRight) + j) + (*p & d->modulo)]++;
            p++;
        }

        // First column initialization
        memset(&H, 0, sizeof(H));
        memset(luc, 0, sizeof(luc));
        if (padLeft)
            histogramMulAdd_sse2(d->radius, &hCoarse[0], H.coarse, bins);
        for (int j = 0; j < (padLeft ? d->radius : 2 * d->radius); j++)
            histogramAdd_sse2(&hCoarse[bins * j], H.coarse, bins);
        for (int k = 0; k < bins; k++)
            histogramMulAdd_sse2(2 * d->radius + 1, &hFine[bins * width * k], &H.fine[k][0], bins);

        for (int j = padLeft ? 0 : d->radius; j < (padRight ? width : width - d->radius); j++) {
            uint16_t sum = 0, * segment;
            int k, b;

            histogramAdd_sse2(&hCoarse[bins * std::min(j + d->radius, width - 1)], H.coarse, bins);

            // Find median at coarse level
            for (k = 0; k < bins; k++) {
                sum += H.coarse[k];
                if (sum > d->t) {
                    sum -= H.coarse[k];
                    break;
                }
            }
            assert(k < 16);

            // Update corresponding histogram segment
            if (luc[k] <= j - d->radius) {
                memset(&H.fine[k], 0, bins * sizeof(uint16_t));
                for (luc[k] = j - d->radius; luc[k] < std::min(j + d->radius + 1, width); luc[k]++)
                    histogramAdd_sse2(&hFine[bins * (width * k + luc[k])], H.fine[k], bins);
                if (luc[k] < j + d->radius + 1) {
                    histogramMulAdd_sse2(j + d->radius + 1 - width, &hFine[bins * (width * k + width - 1)], &H.fine[k][0], bins);
                    luc[k] = j + d->radius + 1;
                }
            } else {
                for (; luc[k] < j + d->radius + 1; luc[k]++) {
                    histogramSub_sse2(&hFine[bins * (width * k + std::max(luc[k] - 2 * d->radius - 1, 0))], H.fine[k], bins);
                    histogramAdd_sse2(&hFine[bins * (width * k + std::min(static_cast<int>(luc[k]), width - 1))], H.fine[k], bins);
                }
            }

            histogramSub_sse2(&hCoarse[bins * std::max(j - d->radius, 0)], H.coarse, bins);

            // Find median in segment
            segment = H.fine[k];
            for (b = 0; b < bins; b++) {
                sum += segment[b];
                if (sum > d->t) {
                    dstp[stride * i + j] = bins * k + b;
                    break;
                }
            }
            assert(b < 16);
        }
    }
}

template void process_sse2<uint8_t, 16>(const uint8_t *, uint8_t *, uint16_t *, uint16_t *, const CTMFData *, const int, const int, const int, const bool, const bool) noexcept;
template void process_sse2<uint16_t, 32>(const uint16_t *, uint16_t *, uint16_t *, uint16_t *, const CTMFData *, const int, const int, const int, const bool, const bool) noexcept;
template void process_sse2<uint16_t, 64>(const uint16_t *, uint16_t *, uint16_t *, uint16_t *, const CTMFData *, const int, const int, const int, const bool, const bool) noexcept;
template void process_sse2<uint16_t, 128>(const uint16_t *, uint16_t *, uint16_t *, uint16_t *, const CTMFData *, const int, const int, const int, const bool, const bool) noexcept;
template void process_sse2<uint16_t, 256>(const uint16_t *, uint16_t *, uint16_t *, uint16_t *, const CTMFData *, const int, const int, const int, const bool, const bool) noexcept;
#endif
