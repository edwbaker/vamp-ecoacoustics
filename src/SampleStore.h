/*
    Ecoacoustic Vamp Plugins

    High-performance implementations of acoustic indices for
    bioacoustics and soundscape ecology analysis.

    (C) Ed Baker 2025. Licensed under GPL (>=2).
*/

// SampleStore.h - the samples of a whole recording, as they were received (floats), in blocks held in chunks.
//
// sh, th and h keep every sample until the end and then take one big FFT. They used to convert each sample to
// double on arrival (8 bytes each, in a std::vector that doubles its buffer when full). Widening a float to a double
// is exact and can be done at the moment the value is used, so the samples are kept as floats (4 bytes each) in a
// RowStore, one row per block, which never moves what it has stored.

#ifndef _SAMPLE_STORE_H_
#define _SAMPLE_STORE_H_

#include "RowStore.h"
#include <algorithm>
#include <cstddef>

class SampleStore
{
public:
    SampleStore() : m_blockSize(1) { }

    // Samples arrive in blocks of blockSize
    void init(size_t blockSize) {
        m_blockSize = blockSize ? blockSize : 1;
        m_rows.init(m_blockSize);
    }

    void clear() { m_rows.clear(); }

    // Storage for the next block
    float* appendBlock() { return m_rows.appendRow(); }

    // Number of samples received (a whole number of blocks)
    size_t size() const { return m_rows.rows() * m_blockSize; }

    float operator[](size_t i) const { return m_rows.row(i / m_blockSize)[i % m_blockSize]; }

    // The number of samples left when the zeros at the end are dropped (the length up to and including the last
    // sample that is not equal to 0; a NaN is not equal to 0)
    size_t sizeWithoutTrailingZeros() const {
        size_t n = size();
        for (size_t row = m_rows.rows(); row > 0 && n > 0; --row) {
            const float* block = m_rows.row(row - 1);
            for (size_t j = m_blockSize; j > 0; --j) {
                if (block[j - 1] != 0.0f) return (row - 1) * m_blockSize + j;
            }
            n = (row - 1) * m_blockSize;
        }
        return 0;
    }

    // The samples first .. first+count-1 in place if they lie in one block, otherwise copied into scratch
    // (which must hold count floats)
    const float* view(size_t first, size_t count, float* scratch) const {
        const size_t row = first / m_blockSize, offset = first % m_blockSize;
        if (offset + count <= m_blockSize) return m_rows.row(row) + offset;
        copyFloats(first, count, scratch);
        return scratch;
    }

    // The samples first .. first+count-1, widened to double
    void copyTo(size_t first, size_t count, double* dst) const {
        size_t row = first / m_blockSize, offset = first % m_blockSize;
        while (count > 0) {
            const size_t n = std::min(count, m_blockSize - offset);
            const float* src = m_rows.row(row) + offset;
            for (size_t i = 0; i < n; ++i) dst[i] = src[i];
            dst += n;
            count -= n;
            ++row;
            offset = 0;
        }
    }

    // Run f(const float* samples, size_t firstIndex, size_t count) over the samples 0 .. n-1, block by block
    template <typename F>
    void forEachBlock(size_t n, F f) const {
        for (size_t row = 0, i = 0; i < n; ++row, i += m_blockSize) {
            f(m_rows.row(row), i, std::min(m_blockSize, n - i));
        }
    }

private:
    void copyFloats(size_t first, size_t count, float* dst) const {
        size_t row = first / m_blockSize, offset = first % m_blockSize;
        while (count > 0) {
            const size_t n = std::min(count, m_blockSize - offset);
            const float* src = m_rows.row(row) + offset;
            std::copy(src, src + n, dst);
            dst += n;
            count -= n;
            ++row;
            offset = 0;
        }
    }

    size_t m_blockSize;
    RowStore m_rows;
};

#endif
