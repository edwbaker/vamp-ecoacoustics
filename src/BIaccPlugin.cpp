/*
    Ecoacoustic Vamp Plugins
    
    High-performance implementations of acoustic indices for 
    bioacoustics and soundscape ecology analysis.
    
    (C) Ed Baker 2025. Licensed under GPL (>=2).
*/

// BIaccPlugin.cpp - Bioacoustic Index (Accumulated)
// Matches soundecology::bioacoustic_index

#include "BIaccPlugin.h"
#include "RealFFT.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace {
// Declared defaults, shared by the constructor and getParameterDescriptors()
const float kDefaultMinFreq = 2000.0f;   // Hz
const float kDefaultMaxFreq = 8000.0f;   // Hz
}

BIaccPlugin::BIaccPlugin(float inputSampleRate) :
    Plugin(inputSampleRate),
    m_channels(0),
    m_blockSize(0),
    m_stepSize(0),
    m_numBins(0),
    m_fft(nullptr),
    m_windowType(Hanning),
    m_minFreq(kDefaultMinFreq),
    m_maxFreq(kDefaultMaxFreq)
{
}

BIaccPlugin::~BIaccPlugin()
{
    if (m_fft) {
        delete m_fft;
        m_fft = nullptr;
    }
}

Vamp::Plugin::InputDomain
BIaccPlugin::getInputDomain() const
{
    return TimeDomain;
}

size_t
BIaccPlugin::getMinChannelCount() const
{
    return 1;
}

size_t
BIaccPlugin::getMaxChannelCount() const
{
    return 2;
}

string
BIaccPlugin::getIdentifier() const
{
    return "bi-acc";
}

string
BIaccPlugin::getName() const
{
    return "Bioacoustic Index (soundecology)";
}

string
BIaccPlugin::getDescription() const
{
    return "Calculates the Bioacoustic Index (BI), matching soundecology implementation.";
}

string
BIaccPlugin::getMaker() const
{
    return "Ecoacoustic-Vamp-Plugins";
}

int
BIaccPlugin::getPluginVersion() const
{
    return 1;
}

string
BIaccPlugin::getCopyright() const
{
    return "(C) Ed Baker 2025. Licensed under GPL (>=2).";
}

Vamp::Plugin::ParameterList
BIaccPlugin::getParameterDescriptors() const
{
    ParameterList list;

    ParameterDescriptor d;
    d.identifier = "minFreq";
    d.name = "Minimum Frequency";
    d.description = "Minimum frequency (Hz)";
    d.unit = "Hz";
    d.minValue = 0;
    d.maxValue = m_inputSampleRate / 2.0f;
    d.defaultValue = kDefaultMinFreq;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "maxFreq";
    d.name = "Maximum Frequency";
    d.description = "Maximum frequency (Hz)";
    d.unit = "Hz";
    d.minValue = 0;
    d.maxValue = m_inputSampleRate / 2.0f;
    d.defaultValue = kDefaultMaxFreq;
    d.isQuantized = false;
    list.push_back(d);

    return list;
}

float
BIaccPlugin::getParameter(string identifier) const
{
    if (identifier == "minFreq") return m_minFreq;
    if (identifier == "maxFreq") return m_maxFreq;
    return 0;
}

void
BIaccPlugin::setParameter(string identifier, float value)
{
    if (identifier == "minFreq") {
        m_minFreq = value;
    } else if (identifier == "maxFreq") {
        m_maxFreq = value;
    }
}

Vamp::Plugin::ProgramList
BIaccPlugin::getPrograms() const
{
    ProgramList list;
    return list;
}

string
BIaccPlugin::getCurrentProgram() const
{
    return "";
}

void
BIaccPlugin::selectProgram(string name)
{
    (void)name;
}

Vamp::Plugin::OutputList
BIaccPlugin::getOutputDescriptors() const
{
    OutputList list;

    OutputDescriptor d;
    d.identifier = "bi";
    d.name = "Bioacoustic Index";
    d.description = "Bioacoustic Index (Left channel or mono)";
    d.unit = "";
    d.hasFixedBinCount = true;
    d.binCount = 1;
    d.hasKnownExtents = false;
    d.isQuantized = false;
    d.sampleType = OutputDescriptor::VariableSampleRate;
    d.hasDuration = false;
    list.push_back(d);

    // Add right channel output for stereo files
    OutputDescriptor d2;
    d2.identifier = "bi_right";
    d2.name = "Bioacoustic Index (Right)";
    d2.description = "Bioacoustic Index (Right channel)";
    d2.unit = "";
    d2.hasFixedBinCount = true;
    d2.binCount = 1;
    d2.hasKnownExtents = false;
    d2.isQuantized = false;
    d2.sampleType = OutputDescriptor::VariableSampleRate;
    d2.hasDuration = false;
    list.push_back(d2);

    return list;
}

size_t
BIaccPlugin::getPreferredBlockSize() const
{
    return 512;
}

size_t
BIaccPlugin::getPreferredStepSize() const
{
    return 512;
}

bool
BIaccPlugin::initialise(size_t channels, size_t stepSize, size_t blockSize)
{
    if (channels < getMinChannelCount() ||
        channels > getMaxChannelCount()) return false;

    m_channels = channels;
    m_blockSize = blockSize;
    m_stepSize = stepSize;
    m_numBins = blockSize / 2;

    // Initialize per-channel vectors
    m_frameCount_ch.assign(m_channels, 0);
    m_globalMax_ch.assign(m_channels, 0.0f);
    m_spectralData_ch.resize(m_channels);
    m_accumulatedSpectrum_ch.resize(m_channels);
    m_pendingPower_ch.assign(m_channels, std::vector<double>(m_numBins, 0.0));
    m_hasPending_ch.assign(m_channels, 0);
    m_pendingPadded_ch.assign(m_channels, 0);
    m_maxBeforeLast_ch.assign(m_channels, 0.0f);

    for (size_t ch = 0; ch < m_channels; ++ch) {
        m_spectralData_ch[ch].clear();
        m_accumulatedSpectrum_ch[ch].assign(m_numBins, 0.0);
    }

    // Initialize FFT
    if (m_fft) delete m_fft;
    m_fft = new InPlaceRealFFT(blockSize);
    
    // Pre-compute window
    m_window.resize(m_blockSize);
    // Symmetric Hanning window to match seewave
    for (size_t i = 0; i < m_blockSize; ++i) {
        m_window[i] = static_cast<float>(0.5 * (1.0 - std::cos(2.0 * M_PI * i / (m_blockSize - 1))));
    }

    // Pre-compute frequency bin range
    float df = (m_inputSampleRate / 2.0f) / static_cast<float>(m_numBins);
    m_rowsWidth = 1.0f / df;
    
    long minBinLong = static_cast<long>(std::floor(m_minFreq * m_rowsWidth)) - 1;
    long maxBinLong = static_cast<long>(std::floor(m_maxFreq * m_rowsWidth));
    
    if (minBinLong < 0) minBinLong = 0;
    
    m_minBin = static_cast<size_t>(minBinLong);
    m_maxBin = static_cast<size_t>(maxBinLong);
    
    if (m_minBin >= m_numBins) m_minBin = m_numBins - 1;
    if (m_maxBin > m_numBins) m_maxBin = m_numBins;

    return true;
}

void
BIaccPlugin::reset()
{
    for (size_t ch = 0; ch < m_channels; ++ch) {
        m_spectralData_ch[ch].clear();
        m_frameCount_ch[ch] = 0;
        m_globalMax_ch[ch] = 0.0f;
        m_accumulatedSpectrum_ch[ch].assign(m_numBins, 0.0);
        m_hasPending_ch[ch] = 0;
        m_pendingPadded_ch[ch] = 0;
        m_maxBeforeLast_ch[ch] = 0.0f;
    }
}

Vamp::Plugin::FeatureSet
BIaccPlugin::process(const float *const *inputBuffers, Vamp::RealTime timestamp)
{
    FeatureSet fs;
    (void)timestamp;
    
    const float* win = m_window.data();
    
    // Process each channel independently: window the block into the FFT's own buffer and transform it straight away
    double* frame = m_fft->buffer();
    for (size_t ch = 0; ch < m_channels; ++ch) {
        const float* src = inputBuffers[ch];
        for (size_t i = 0; i < m_blockSize; ++i) {
            frame[i] = static_cast<double>(src[i]) * win[i];
        }

        // A block ending in zeros is the host's padding of an incomplete last window (the rule ACIaccPlugin uses)
        processFrame(ch, src[m_blockSize - 1] == 0.0f);
    }
    
    return fs;
}

void BIaccPlugin::processFrame(size_t channel, bool padded)
{
    double* accSpectrum = m_accumulatedSpectrum_ch[channel].data();
    double* pending = m_pendingPower_ch[channel].data();
    float& globalMax = m_globalMax_ch[channel];

    // The previous frame was not the last one: add it to the sums now (the additions are in the same order as before)
    if (m_hasPending_ch[channel]) {
        for (size_t b = 0; b < m_numBins; ++b) {
            accSpectrum[b] += pending[b];
        }
        m_frameCount_ch[channel]++;
    }
    m_maxBeforeLast_ch[channel] = globalMax;

    // Perform FFT
    m_fft->forward();

    // Compute magnitudes; the powers wait in pending until the next frame or getRemainingFeatures()
    // Bin 0 (DC)
    double real0 = m_fft->realOfBin0();
    double mag0 = std::abs(real0);
    if (static_cast<float>(mag0) > globalMax) globalMax = static_cast<float>(mag0);
    pending[0] = mag0 * mag0; // Power

    // Bins 1 to N/2 - 1. The global maximum is float(sqrt(magSq)) over all bins; float(sqrt(x)) never decreases as
    // x grows, so it is the same as float(sqrt(largest magSq)), which needs one sqrt per frame instead of one per bin.
    double maxSq = 0.0;
    m_fft->forBins(1, m_numBins, [pending, &maxSq](size_t i, double real, double imag) {
        double magSq = real * real + imag * imag;

        if (magSq > maxSq) maxSq = magSq;

        pending[i] = magSq;
    });
    if (static_cast<float>(std::sqrt(maxSq)) > globalMax) {
        globalMax = static_cast<float>(std::sqrt(maxSq));
    }

    m_hasPending_ch[channel] = 1;
    m_pendingPadded_ch[channel] = padded;
}

Vamp::Plugin::FeatureSet
BIaccPlugin::getRemainingFeatures()
{
    FeatureSet fs;
    
    // Calculate BI for each channel
    for (size_t ch = 0; ch < m_channels; ++ch) {
        // seewave::spectro uses only complete windows: the newest frame counts unless the host padded its block with
        // zeros, or it is the only frame
        if (m_hasPending_ch[ch]) {
            if (!m_pendingPadded_ch[ch] || m_frameCount_ch[ch] == 0) {
                const double* pending = m_pendingPower_ch[ch].data();
                for (size_t b = 0; b < m_numBins; ++b) {
                    m_accumulatedSpectrum_ch[ch][b] += pending[b];
                }
                m_frameCount_ch[ch]++;
            } else {
                m_globalMax_ch[ch] = m_maxBeforeLast_ch[ch];
            }
            m_hasPending_ch[ch] = 0;
        }

        if (m_frameCount_ch[ch] == 0) {
            Feature f;
            f.hasTimestamp = true;
            f.timestamp = Vamp::RealTime::frame2RealTime(0, static_cast<unsigned int>(m_inputSampleRate));
            f.values.push_back(0.0f);
            fs[static_cast<int>(ch)].push_back(f);
            continue;
        }

        // Calculate Mean Spectrum (Power)
        std::vector<double> meanSpectrum(m_numBins);
        
        for (size_t b = 0; b < m_numBins; ++b) {
            meanSpectrum[b] = m_accumulatedSpectrum_ch[ch][b] / m_frameCount_ch[ch];
        }
        
        // Convert to dB Normalized (Max 0)
        float maxAmp = m_globalMax_ch[ch];
        if (maxAmp <= 1e-9f) maxAmp = 1.0f; 
        double maxPower = maxAmp * maxAmp;

        std::vector<double> meanSpectrumDB(m_numBins);
        for (size_t b = 0; b < m_numBins; ++b) {
            double val = meanSpectrum[b] / maxPower;
            // No floor, as in soundecology. Only an exactly zero power, where log10 gives -inf (and soundecology a
            // NaN), is held at 1e-300 (-3000 dB), so that silence still gives a BI of 0
            if (val <= 0.0) val = 1e-300;
            meanSpectrumDB[b] = 10.0 * std::log10(val);
        }

        // Extract Segment and Normalize
        size_t minBin = m_minBin;
        size_t maxBin = m_maxBin;
        
        if (minBin >= maxBin) {
            Feature f;
            f.hasTimestamp = true;
            f.timestamp = Vamp::RealTime::frame2RealTime(0, static_cast<unsigned int>(m_inputSampleRate));
            f.values.push_back(0.0f);
            fs[static_cast<int>(ch)].push_back(f);
            continue;
        }

        double minVal = 1e9;
        for (size_t b = minBin; b < maxBin; ++b) {
            if (meanSpectrumDB[b] < minVal) minVal = meanSpectrumDB[b];
        }
        
        // Calculate Area
        double bi = 0.0;
        double rowsWidth = static_cast<double>(m_rowsWidth);
        for (size_t b = minBin; b < maxBin; ++b) {
            bi += (meanSpectrumDB[b] - minVal) * rowsWidth;
        }
        
        Feature f;
        f.hasTimestamp = true;
        f.timestamp = Vamp::RealTime::frame2RealTime(0, static_cast<unsigned int>(m_inputSampleRate));
        f.values.push_back(static_cast<float>(bi));
        
        fs[static_cast<int>(ch)].push_back(f);
    }
    
    return fs;
}
