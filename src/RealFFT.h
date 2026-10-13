#ifndef ECOACOUSTICS_REALFFT_H
#define ECOACOUSTICS_REALFFT_H

#include <algorithm>
#include <cstddef>
#include <vector>

// A forward real FFT that works on the buffer it owns. The caller writes the (windowed) samples straight into
// buffer(), calls forward(), and reads the bins back with forBins(). This is the computation of
// Vamp::FFTReal::forward -- the same pocketfft plan, run in place on the same n values -- without its copy of
// the input into a work array and its unpacking of the result into an interleaved complex array. Every value
// read back is the same double that FFTReal::forward would have written, so results are bit-identical.
// (The transform itself is compiled in FFT.cpp, like FFTReal's, so that the plugins do not each carry a copy.)
class InPlaceRealFFT
{
public:
    explicit InPlaceRealFFT(size_t n);
    ~InPlaceRealFFT();
    InPlaceRealFFT(const InPlaceRealFFT&) = delete;
    InPlaceRealFFT& operator=(const InPlaceRealFFT&) = delete;

    size_t size() const { return m_n; }

    // n doubles to fill with the samples before forward(); holds pocketfft's half-complex result afterwards
    double* buffer() { return m_buf.data(); }

    void forward();

    // The real part of bin 0 (the sum of the samples) after forward()
    double realOfBin0() const { return m_buf[0]; }

    // Call f(k, real, imag) for the bins first <= k < last of the last forward(), in increasing k, with the
    // values FFTReal::forward returned for them: bin 0 has imaginary part +0.0, so does the Nyquist bin of an
    // even length, and a bin above n/2 reads as (0.0, 0.0).
    template <class F>
    void forBins(size_t first, size_t last, F&& f) const {
        const double* h = m_buf.data();
        const size_t n = m_n;
        size_t k = first;
        if (k == 0 && last > 0) {
            f(k, h[0], 0.0);
            k = 1;
        }
        // bins 1 .. (n-1)/2 have a real and an imaginary part: h[2k-1], h[2k]
        const size_t interiorEnd = std::min(last, (n > 0 ? (n - 1) / 2 : 0) + 1);
        for (; k < interiorEnd; ++k) {
            f(k, h[2 * k - 1], h[2 * k]);
        }
        if (k < last && n >= 2 && n % 2 == 0 && k == n / 2) {
            f(k, h[n - 1], 0.0);
            ++k;
        }
        for (; k < last; ++k) {
            f(k, 0.0, 0.0);
        }
    }

private:
    struct Plan;
    size_t m_n;
    std::vector<double> m_buf;
    Plan* m_plan;
};

#endif
