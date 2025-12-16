// THPlugin.cpp - Temporal Entropy (Ht)
// Matches seewave::th(env(wave)) behavior
// Uses PocketFFT for Hilbert transform on arbitrary sizes

#include "THPlugin.h"
#include "../ext/pocketfft_hdronly.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <complex>

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
    m_inputBuffer.clear();
    m_inputBuffer.reserve(blockSize * 1000);

    return true;
}

void THPlugin::reset()
{
    m_inputBuffer.clear();
    m_totalSamples = 0;
}

THPlugin::FeatureSet THPlugin::process(const float *const *inputBuffers, Vamp::RealTime timestamp)
{
    (void)timestamp;
    if (m_stepSize == 0 || m_blockSize == 0) return FeatureSet();

    // Accumulate ALL samples - we'll do Hilbert transform on the full signal at the end
    size_t oldSize = m_inputBuffer.size();
    m_inputBuffer.resize(oldSize + m_blockSize);
    double* dst = m_inputBuffer.data() + oldSize;
    const float* src = inputBuffers[0];
    for (size_t i = 0; i < m_blockSize; ++i) {
        dst[i] = static_cast<double>(src[i]);
    }

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
    
    size_t n = m_inputBuffer.size();
    if (n == 0) return fs;

    // Remove trailing zeros that may be padding from the last block
    while (n > 0 && m_inputBuffer[n-1] == 0.0) {
        n--;
    }
    // Add back a small buffer in case file actually ends with zeros
    size_t removed = m_inputBuffer.size() - n;
    if (removed > 0 && removed < m_blockSize) {
        m_inputBuffer.resize(n);
    } else {
        n = m_inputBuffer.size();
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
    
    // Pad with zeros (centered, like seewave)
    std::vector<double> padded(n_padded, 0.0);
    size_t nzp = n_padded - n;
    size_t pad_start = nzp / 2;
    for (size_t i = 0; i < n; ++i) {
        padded[pad_start + i] = m_inputBuffer[i];
    }
    
    // Forward FFT using PocketFFT (real to complex)
    std::vector<std::complex<double>> fftOut(n_padded / 2 + 1);
    
    pocketfft::shape_t shape = {n_padded};
    pocketfft::stride_t stride_in = {sizeof(double)};
    pocketfft::stride_t stride_out = {sizeof(std::complex<double>)};
    pocketfft::shape_t axes = {0};
    
    pocketfft::r2c(shape, stride_in, stride_out, axes, 
                   pocketfft::FORWARD, padded.data(), fftOut.data(), 1.0);
    
    std::vector<std::complex<double>> analyticFFT(n_padded);
    
    // Index 0 (DC)
    analyticFFT[0] = fftOut[0];
    
    // Positive frequencies (1 to n/2-1) - doubled
    for (size_t k = 1; k < n_padded / 2; ++k) {
        analyticFFT[k] = 2.0 * fftOut[k];
    }
    
    // Nyquist (n/2) - stays as is
    analyticFFT[n_padded / 2] = fftOut[n_padded / 2];
    
    // Negative frequencies (n/2+1 to n-1) - zero
    for (size_t k = n_padded / 2 + 1; k < n_padded; ++k) {
        analyticFFT[k] = std::complex<double>(0.0, 0.0);
    }
    
    // Inverse FFT to get analytic signal
    std::vector<std::complex<double>> analytic(n_padded);
    
    pocketfft::stride_t stride_c = {sizeof(std::complex<double>)};
    pocketfft::c2c(shape, stride_c, stride_c, axes,
                   pocketfft::BACKWARD, analyticFFT.data(), analytic.data(), 1.0 / n_padded);
    
    // Extract envelope (magnitude) and trim back to original length
    std::vector<double> envelope(n_original);
    size_t trim_start = nzp / 2;
    for (size_t i = 0; i < n_original; ++i) {
        envelope[i] = std::abs(analytic[trim_start + i]);
    }
    
    // Now compute TH on the envelope
    size_t N = n_original;
    const double epsilon = 1e-7;
    
    // Compute sum for PMF
    double sumEnv = 0.0;
    for (size_t i = 0; i < N; ++i) {
        sumEnv += envelope[i];
    }
    
    if (sumEnv == 0.0) {
        Feature f;
        f.hasTimestamp = true;
        f.timestamp = Vamp::RealTime::zeroTime;
        f.values.push_back(0.0f);
        fs[0].push_back(f);
        return fs;
    }
    
    // Create PMF and replace zeros with epsilon
    std::vector<double> pmf(N);
    for (size_t i = 0; i < N; ++i) {
        pmf[i] = envelope[i] / sumEnv;
        if (pmf[i] == 0.0) pmf[i] = epsilon;
    }
    
    // Compute Shannon entropy
    double entropySum = 0.0;
    for (size_t i = 0; i < N; ++i) {
        entropySum += pmf[i] * std::log(pmf[i]);
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
