/*
    Ecoacoustic Vamp Plugins
    
    High-performance implementations of acoustic indices for 
    bioacoustics and soundscape ecology analysis.
    
    (C) Ed Baker 2025. Licensed under GPL (>=2).
*/

// SHPlugin.cpp - Spectral Entropy
// Matches seewave::sh(spec(wave)) behavior
// Uses PocketFFT for efficient FFT on arbitrary sizes

#include "SHPlugin.h"
#include "../ext/pocketfft_hdronly.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <complex>
#include <memory>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

SHPlugin::SHPlugin(float inputSampleRate) :
    Plugin(inputSampleRate),
    m_channels(0),
    m_blockSize(0),
    m_stepSize(0),
    m_numBins(0),
    m_frameCount(0),
    m_globalMax(0.0f),
    m_fft(nullptr),
    m_batchSize(256),
    m_windowType(Hanning),
    m_minFreq(0.0f),
    m_maxFreq(0.0f)
{
}

SHPlugin::~SHPlugin()
{
    if (m_fft) {
        delete m_fft;
        m_fft = nullptr;
    }
}

Vamp::Plugin::InputDomain
SHPlugin::getInputDomain() const
{
    return TimeDomain;
}

size_t
SHPlugin::getMinChannelCount() const
{
    return 1;
}

size_t
SHPlugin::getMaxChannelCount() const
{
    return 2;
}

string
SHPlugin::getIdentifier() const
{
    return "sh";
}

string
SHPlugin::getName() const
{
    return "Spectral Entropy (seewave)";
}

string
SHPlugin::getDescription() const
{
    return "Calculates the Spectral Entropy (SH) of a signal, matching seewave implementation.";
}

string
SHPlugin::getMaker() const
{
    return "Ecoacoustic-Vamp-Plugins";
}

int
SHPlugin::getPluginVersion() const
{
    return 1;
}

string
SHPlugin::getCopyright() const
{
    return "(C) Ed Baker 2025. Licensed under GPL (>=2).";
}

Vamp::Plugin::ParameterList
SHPlugin::getParameterDescriptors() const
{
    ParameterList list;
    return list;
}

float
SHPlugin::getParameter(string identifier) const
{
    (void)identifier;
    return 0;
}

void
SHPlugin::setParameter(string identifier, float value)
{
    (void)identifier;
    (void)value;
}

Vamp::Plugin::ProgramList
SHPlugin::getPrograms() const
{
    ProgramList list;
    return list;
}

string
SHPlugin::getCurrentProgram() const
{
    return "";
}

void
SHPlugin::selectProgram(string name)
{
    (void)name;
}

Vamp::Plugin::OutputList
SHPlugin::getOutputDescriptors() const
{
    OutputList list;

    OutputDescriptor d;
    d.identifier = "sh";
    d.name = "Spectral Entropy";
    d.description = "Spectral Entropy (Shannon entropy of the spectral distribution)";
    d.unit = "";
    d.hasFixedBinCount = true;
    d.binCount = 1;
    d.hasKnownExtents = false;
    d.isQuantized = false;
    d.sampleType = OutputDescriptor::VariableSampleRate;
    d.hasDuration = false;
    list.push_back(d);

    return list;
}

size_t
SHPlugin::getPreferredBlockSize() const
{
    return 512;
}

size_t
SHPlugin::getPreferredStepSize() const
{
    return 512;
}

bool
SHPlugin::initialise(size_t channels, size_t stepSize, size_t blockSize)
{
    if (channels < getMinChannelCount() ||
        channels > getMaxChannelCount()) return false;

    m_channels = channels;
    m_blockSize = blockSize;
    m_stepSize = stepSize;
    m_frameCount = 0;
    m_globalMax = 0.0f;
    m_spectralData.clear();
    m_numBins = 0;

    // Accumulate all samples (the one FFT is made at the end, directly with pocketfft)
    m_samples.init(m_blockSize);
    
    // Hanning window for STFT frames
    m_window.resize(m_blockSize);
    for (size_t i = 0; i < m_blockSize; ++i) {
        m_window[i] = 0.5 * (1.0 - std::cos(2.0 * M_PI * i / (m_blockSize - 1)));
    }
    
    m_accumulatedSpectrum.clear();
    m_totalSamples = 0;

    return true;
}

void
SHPlugin::reset()
{
    m_spectralData.clear();
    m_samples.clear();
    m_frameCount = 0;
    m_globalMax = 0.0f;
    m_accumulatedSpectrum.clear();
    m_totalSamples = 0;
}

Vamp::Plugin::FeatureSet
SHPlugin::process(const float *const *inputBuffers, Vamp::RealTime timestamp)
{
    FeatureSet fs;
    
    // Track actual sample position from timestamp
    // This gives us the real sample count without block padding
    size_t frameNumber = Vamp::RealTime::realTime2Frame(timestamp, static_cast<unsigned int>(m_inputSampleRate));
    m_totalSamples = frameNumber + m_blockSize;
    
    // Accumulate ALL samples - we'll do one big FFT at the end
    std::copy(inputBuffers[0], inputBuffers[0] + m_blockSize, m_samples.appendBlock());
    
    return fs;
}

void SHPlugin::processBatch(size_t numFrames)
{
    // Not used
    (void)numFrames;
}

Vamp::Plugin::FeatureSet
SHPlugin::getRemainingFeatures()
{
    FeatureSet fs;
    
    const size_t received = m_samples.size();
    size_t n = received;
    if (n == 0) return fs;

    // Remove trailing zeros that may be padding from the last block
    // The host may pad the last block with zeros to fill it
    n = m_samples.sizeWithoutTrailingZeros();
    // Add back a small buffer in case file actually ends with zeros
    // This is a heuristic - check if we removed more than a block worth
    size_t removed = received - n;
    if (removed > 0 && removed < m_blockSize) {
        // Likely actual padding, keep the trimmed size
    } else {
        // Either no zeros removed or too many - keep original
        n = received;
    }

    // Apply Hanning window to entire signal (seewave default), widening each sample to double as it is used
    std::unique_ptr<double[]> data(new double[n]);
    m_samples.forEachBlock(n, [&](const float* src, size_t first, size_t count) {
        for (size_t j = 0; j < count; ++j) {
            const size_t i = first + j;
            double w = 0.5 * (1.0 - std::cos(2.0 * M_PI * i / (n - 1)));
            data[i] = static_cast<double>(src[j]) * w;
        }
    });
    m_samples.clear();                              // no longer needed

    // Perform the real-to-complex FFT with the PocketFFT plan that pocketfft::r2c uses, on the array itself.
    // r2c copies its input line into a work array of its own, runs this same exec on that and unpacks the
    // result; running it here, in place, gives the same numbers without the second array. The result is in
    // PocketFFT's half-complex order: r0, r1, i1, r2, i2, ... (with the real-only Nyquist value last for even n).
    {
        pocketfft::detail::pocketfft_r<double> plan(n);
        plan.exec(data.get(), 1.0, true);
    }

    // Compute magnitude spectrum: bin 0 is (r0, 0), bin i is (r_i, i_i) = (data[2i-1], data[2i])
    m_numBins = n / 2;
    m_accumulatedSpectrum.resize(m_numBins);

    const double* spec = data.get();
    if (m_numBins > 0) {
        m_accumulatedSpectrum[0] = 2.0 * std::abs(std::complex<double>(spec[0], 0.0));
    }
    for (size_t i = 1; i < m_numBins; ++i) {
        m_accumulatedSpectrum[i] = 2.0 * std::abs(std::complex<double>(spec[2 * i - 1], spec[2 * i]));
    }
    data.reset();                                   // no longer needed

    m_frameCount = 1;

    double sh = computeSH();

    // The spectrum is not needed any more: release it now rather than with the plugin
    std::vector<double>().swap(m_accumulatedSpectrum);
    m_numBins = 0;
    
    Feature f;
    f.hasTimestamp = true;
    f.timestamp = Vamp::RealTime::frame2RealTime(0, static_cast<unsigned int>(m_inputSampleRate));
    f.values.push_back(static_cast<float>(sh));
    
    fs[0].push_back(f);
    
    return fs;
}

double SHPlugin::computeSH() const
{
    if (m_numBins == 0) return 0.0;

    size_t N = m_numBins;
    const double epsilon = 1e-7;
    
    // First normalize to max=1 and replace zeros with epsilon
    double maxSpectrum = 0.0;
    for (size_t i = 0; i < N; ++i) {
        if (m_accumulatedSpectrum[i] > maxSpectrum) {
            maxSpectrum = m_accumulatedSpectrum[i];
        }
    }
    if (maxSpectrum == 0.0) return 0.0;
    
    // Create normalized spectrum with zeros replaced by epsilon
    std::vector<double> normalized(N);
    double sumNormalized = 0.0;
    for (size_t i = 0; i < N; ++i) {
        double val = m_accumulatedSpectrum[i] / maxSpectrum;
        if (val == 0.0) val = epsilon;
        normalized[i] = val;
        sumNormalized += val;
    }
    
    if (sumNormalized == 0.0) return 0.0;
    
    // Compute Shannon entropy on PMF
    double entropySum = 0.0;
    for (size_t i = 0; i < N; ++i) {
        double p = normalized[i] / sumNormalized; // PMF
        entropySum += p * std::log(p); // Natural log like seewave
    }
    
    // Normalize by log(N)
    double sh = -entropySum / std::log(static_cast<double>(N));
    
    return sh;
}

