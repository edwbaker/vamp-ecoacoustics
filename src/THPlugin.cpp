/*
    Ecoacoustic Vamp Plugins
    
    High-performance implementations of acoustic indices for 
    bioacoustics and soundscape ecology analysis.
    
    (C) Ed Baker 2025. Licensed under GPL (>=2).
*/

// THPlugin.cpp - Temporal Entropy (Ht)
// Matches seewave::th(env(wave)) behavior
// Uses PocketFFT for Hilbert transform on arbitrary sizes

#include "THPlugin.h"
#include "../ext/pocketfft_hdronly.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <complex>
#include <memory>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using std::string;
using std::vector;

THPlugin::THPlugin(float inputSampleRate) :
    Plugin(inputSampleRate),
    m_blockSize(0),
    m_stepSize(0),
    m_totalSamples(0)
{
}

THPlugin::~THPlugin()
{
}

string
THPlugin::getIdentifier() const
{
    return "th";
}

string
THPlugin::getName() const
{
    return "Temporal Entropy (Ht)";
}

string
THPlugin::getDescription() const
{
    return "Calculates the Temporal Entropy (Ht) matching seewave::th(env()).";
}

string
THPlugin::getMaker() const
{
    return "Ecoacoustic-Vamp-Plugins";
}

int
THPlugin::getPluginVersion() const
{
    return 1;
}

string
THPlugin::getCopyright() const
{
    return "(C) Ed Baker 2025. Licensed under GPL (>=2).";
}

THPlugin::InputDomain
THPlugin::getInputDomain() const
{
    return TimeDomain;
}

size_t
THPlugin::getPreferredBlockSize() const
{
    return 512;
}

size_t
THPlugin::getPreferredStepSize() const
{
    return 512;
}

size_t
THPlugin::getMinChannelCount() const
{
    return 1;
}

size_t
THPlugin::getMaxChannelCount() const
{
    return 2;
}

THPlugin::ParameterList THPlugin::getParameterDescriptors() const {
    return ParameterList();
}

float THPlugin::getParameter(string identifier) const {
    (void)identifier;
    return 0;
}

void THPlugin::setParameter(string identifier, float value) {
    (void)identifier;
    (void)value;
}

THPlugin::ProgramList THPlugin::getPrograms() const { return ProgramList(); }
string THPlugin::getCurrentProgram() const { return ""; }
void THPlugin::selectProgram(string name) { (void)name; }

THPlugin::OutputList THPlugin::getOutputDescriptors() const
{
    OutputList list;

    OutputDescriptor d;
    d.identifier = "th";
    d.name = "Temporal Entropy (Ht)";
    d.description = "Temporal entropy of the amplitude envelope";
    d.unit = "";
    d.hasFixedBinCount = true;
    d.binCount = 1;
    d.hasKnownExtents = false;
    d.isQuantized = false;
    d.sampleType = OutputDescriptor::VariableSampleRate; // Summary feature
    d.hasDuration = false;
    list.push_back(d);

    return list;
}

bool THPlugin::initialise(size_t channels, size_t stepSize, size_t blockSize)
{
    if (channels < getMinChannelCount() || channels > getMaxChannelCount()) return false;

    m_blockSize = blockSize;
    m_stepSize = stepSize;
    m_totalSamples = 0;
    
    // Accumulate all samples for full-signal processing
    m_samples.init(blockSize);

    return true;
}

void THPlugin::reset()
{
    m_samples.clear();
    m_totalSamples = 0;
}

THPlugin::FeatureSet THPlugin::process(const float *const *inputBuffers, Vamp::RealTime timestamp)
{
    (void)timestamp;
    if (m_stepSize == 0 || m_blockSize == 0) return FeatureSet();

    // Accumulate ALL samples - we'll do Hilbert transform on the full signal at the end
    std::copy(inputBuffers[0], inputBuffers[0] + m_blockSize, m_samples.appendBlock());

    return FeatureSet();
}

double THPlugin::computeTH() const
{
    // This is now computed inline in getRemainingFeatures
    return 0.0;
}

THPlugin::FeatureSet THPlugin::getRemainingFeatures()
{
    FeatureSet fs;
    
    const size_t received = m_samples.size();
    size_t n = received;
    if (n == 0) return fs;

    // Remove trailing zeros that may be padding from the last block
    n = m_samples.sizeWithoutTrailingZeros();
    // Add back a small buffer in case file actually ends with zeros
    size_t removed = received - n;
    if (removed > 0 && removed < m_blockSize) {
        // keep the trimmed size
    } else {
        n = received;
    }
    
    if (n == 0) return fs;

    // seewave::hilbert pads to power of 2, does FFT, then trims back
    // We need to do the same to match exactly
    
    size_t n_original = n;
    
    // Check if already power of 2
    double log_n = std::log2(static_cast<double>(n));
    size_t n_padded = n;
    if (log_n != std::floor(log_n)) {
        size_t p = static_cast<size_t>(std::ceil(log_n));
        n_padded = static_cast<size_t>(1) << p;
    }
    
    // Pad with zeros (centered, like seewave). Only the margins need zeroing: the samples overwrite the rest.
    size_t nzp = n_padded - n;
    size_t pad_start = nzp / 2;
    std::unique_ptr<double[]> padded(new double[n_padded]);
    std::fill(padded.get(), padded.get() + pad_start, 0.0);
    m_samples.copyTo(0, n, padded.get() + pad_start);
    std::fill(padded.get() + pad_start + n, padded.get() + n_padded, 0.0);
    // The accumulated samples are not needed again: release them before the large FFT buffers are allocated.
    m_samples.clear();

    // Forward FFT using the PocketFFT plan that pocketfft::r2c uses, run in place on the padded signal (r2c would
    // copy it into a work array of its own first, then unpack the result: same exec, same numbers).
    // The result is in half-complex order: r0, r1, i1, r2, i2, ... (real-only Nyquist value last for even length).
    {
        pocketfft::detail::pocketfft_r<double> plan(n_padded);
        plan.exec(padded.get(), 1.0, true);
    }

    // The analytic spectrum: complex values for the bins 0..n/2 from the half-complex result, with the positive
    // frequencies (1 to n/2-1) doubled, DC and Nyquist as they are and the negative frequencies zero.
    // (2.0 * complex is 2.0 * real and 2.0 * imaginary, as here.) Every element is written once, so the array is
    // not zero-filled first: it is raw memory from operator new, in which the std::complex<double> objects (an
    // implicit-lifetime type) exist, and their real and imaginary parts are written through the double view that
    // std::complex guarantees (reinterpret_cast<double(&)[2]> of an element), not as a separate array of doubles.
    struct RawMemory { void operator()(void* p) const { ::operator delete(p); } };
    std::unique_ptr<std::complex<double>, RawMemory> analytic(
        static_cast<std::complex<double>*>(::operator new(n_padded * sizeof(std::complex<double>))));
    std::complex<double>* analyticFFT = analytic.get();
    {
        const double* w = padded.get();
        double* a = reinterpret_cast<double*>(analyticFFT);
        const size_t half = n_padded / 2;

        // Index 0 (DC) stays as it is.
        a[0] = w[0];
        a[1] = 0.0;

        // Positive frequencies (1 to n/2-1) - doubled
        for (size_t k = 1; k < half; ++k) {
            a[2 * k] = 2.0 * w[2 * k - 1];
            a[2 * k + 1] = 2.0 * w[2 * k];
        }

        // Nyquist (n/2) stays as it is.
        if (half >= 1) {
            if (n_padded % 2 == 0) {
                a[2 * half] = w[n_padded - 1];
                a[2 * half + 1] = 0.0;
            } else {
                a[2 * half] = w[n_padded - 2];
                a[2 * half + 1] = w[n_padded - 1];
            }
        }

        // Negative frequencies (n/2+1 to n-1) - zero
        for (size_t k = half + 1; k < n_padded; ++k) {
            a[2 * k] = 0.0;
            a[2 * k + 1] = 0.0;
        }
    }
    padded.reset();

    // Inverse FFT, in place, to get the analytic signal
    pocketfft::shape_t shape = {n_padded};
    pocketfft::shape_t axes = {0};
    pocketfft::stride_t stride_c = {sizeof(std::complex<double>)};
    pocketfft::c2c(shape, stride_c, stride_c, axes,
                   pocketfft::BACKWARD, analyticFFT, analyticFFT, 1.0 / n_padded);

    // Extract envelope (magnitude) and trim back to original length. The envelope is written over the start of
    // the analytic array: the value for index i goes to byte offset 8*i, which belongs to element i/2, and that
    // has been read already (or is below trim_start and not used at all) when index i is reached.
    // The sum for the PMF is accumulated in the same loop, in the same order.
    double* envelope = reinterpret_cast<double*>(analyticFFT);
    size_t trim_start = nzp / 2;
    double sumEnv = 0.0;
    for (size_t i = 0; i < n_original; ++i) {
        double e = std::abs(analyticFFT[trim_start + i]);
        envelope[i] = e;
        sumEnv += e;
    }

    // Now compute TH on the envelope
    size_t N = n_original;
    const double epsilon = 1e-7;

    if (sumEnv == 0.0) {
        Feature f;
        f.hasTimestamp = true;
        f.timestamp = Vamp::RealTime::zeroTime;
        f.values.push_back(0.0f);
        fs[0].push_back(f);
        return fs;
    }
    
    // PMF (zeros replaced with epsilon) and Shannon entropy, accumulated in the same order as before
    double entropySum = 0.0;
    for (size_t i = 0; i < N; ++i) {
        double p = envelope[i] / sumEnv;
        if (p == 0.0) p = epsilon;
        entropySum += p * std::log(p);
    }
    
    double Ht = -entropySum / std::log(static_cast<double>(N));

    // Output Ht
    Feature f;
    f.hasTimestamp = true;
    f.timestamp = Vamp::RealTime::zeroTime;
    f.values.push_back(static_cast<float>(Ht));
    fs[0].push_back(f);

    return fs;
}
