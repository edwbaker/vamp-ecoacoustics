/*
    Ecoacoustic Vamp Plugins

    High-performance implementations of acoustic indices for
    bioacoustics and soundscape ecology analysis.

    (C) Ed Baker 2025. Licensed under GPL (>=2).
*/

// RowStore.h - append-only store of equal-length rows, kept in fixed-size chunks.
//
// The plugins that keep every frame of a recording until getRemainingFeatures() used to append to one
// std::vector. That doubles its buffer when it is full (copying everything stored so far, with the old and the
// new buffer alive together) and zero-fills every element before it is overwritten. Here a row, once written, is
// never moved, and new storage is not initialised. Reading a row back gives exactly what was written.

#ifndef _ROW_STORE_H_
#define _ROW_STORE_H_

#include <cstddef>
#include <memory>
#include <vector>

template <typename T>
class BasicRowStore
{
public:
    BasicRowStore() : m_rowLen(1), m_rowsPerChunk(1), m_rows(0) { }

    // Rows of rowLen values; chunks of just under chunkBytes (the slack leaves room for the allocator's own
    // header, so that a chunk of a power-of-two size does not spill into one more page that is then touched)
    void init(size_t rowLen, size_t chunkBytes = 1048576) {
        m_chunks.clear();
        m_rows = 0;
        m_rowLen = rowLen ? rowLen : 1;
        m_rowsPerChunk = (chunkBytes - 256) / (m_rowLen * sizeof(T));
        if (m_rowsPerChunk < 1) m_rowsPerChunk = 1;
    }

    void clear() { m_chunks.clear(); m_rows = 0; }

    size_t rows() const { return m_rows; }
    size_t rowLength() const { return m_rowLen; }
    // Rows are contiguous within a chunk: row(i) for i a multiple of rowsPerChunk() starts a run of
    // min(rowsPerChunk(), rows() - i) * rowLength() values
    size_t rowsPerChunk() const { return m_rowsPerChunk; }

    // Storage for the next row (rowLength() values, contents not initialised)
    T* appendRow() {
        const size_t chunk = m_rows / m_rowsPerChunk;
        if (chunk == m_chunks.size()) {
            m_chunks.emplace_back(new T[m_rowLen * m_rowsPerChunk]);
        }
        T* r = m_chunks[chunk].get() + (m_rows % m_rowsPerChunk) * m_rowLen;
        ++m_rows;
        return r;
    }

    T* row(size_t i) {
        return m_chunks[i / m_rowsPerChunk].get() + (i % m_rowsPerChunk) * m_rowLen;
    }
    const T* row(size_t i) const {
        return m_chunks[i / m_rowsPerChunk].get() + (i % m_rowsPerChunk) * m_rowLen;
    }

private:
    size_t m_rowLen;
    size_t m_rowsPerChunk;
    size_t m_rows;
    std::vector<std::unique_ptr<T[]>> m_chunks;
};

typedef BasicRowStore<float> RowStore;
typedef BasicRowStore<double> RowStoreD;

#endif
