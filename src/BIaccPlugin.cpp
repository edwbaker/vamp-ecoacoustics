// BIaccPlugin.cpp - Bioacoustic Index (Accumulated)
// Matches soundecology::bioacoustic_index

#include "BIaccPlugin.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

BIaccPlugin::BIaccPlugin(float inputSampleRate) :
    Plugin(inputSampleRate),
    m_channels(0),
    m_blockSize(0),
    m_stepSize(0),
    m_numBins(0),
    m_fft(nullptr),
    m_batchSize(256),
    m_windowType(Hanning),
    m_minFreq(2000.0f),
    m_maxFreq(8000.0f)
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
    d.defaultValue = 2000.0f;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "maxFreq";
    d.name = "Maximum Frequency";
    d.description = "Maximum frequency (Hz)";
    d.unit = "Hz";
    d.minValue = 0;
    d.maxValue = m_inputSampleRate / 2.0f;
    d.defaultValue = 8000.0f;
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
    m_inputBuffer_ch.resize(m_channels);
    m_spectralData_ch.resize(m_channels);
    m_accumulatedSpectrum_ch.resize(m_channels);
    
    for (size_t ch = 0; ch < m_channels; ++ch) {
        m_inputBuffer_ch[ch].reserve(blockSize * m_batchSize);
        m_spectralData_ch[ch].clear();
        m_accumulatedSpectrum_ch[ch].assign(m_numBins, 0.0);
    }

    // Initialize FFT
    if (m_fft) delete m_fft;
    m_fft = new Vamp::FFTReal(static_cast<unsigned int>(blockSize));
    m_fftOut.resize(blockSize + 2);
    
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
        m_inputBuffer_ch[ch].clear();
        m_frameCount_ch[ch] = 0;
        m_globalMax_ch[ch] = 0.0f;
        m_accumulatedSpectrum_ch[ch].assign(m_numBins, 0.0);
    }
}

Vamp::Plugin::FeatureSet
BIaccPlugin::process(const float *const *inputBuffers, Vamp::RealTime timestamp)
{
    FeatureSet fs;
    (void)timestamp;
    
    const float* win = m_window.data();
    
    // Process each channel independently
    for (size_t ch = 0; ch < m_channels; ++ch) {
        // Apply window and append to input buffer for this channel
        size_t oldSize = m_inputBuffer_ch[ch].size();
        m_inputBuffer_ch[ch].resize(oldSize + m_blockSize);
        double* dst = m_inputBuffer_ch[ch].data() + oldSize;
        const float* src = inputBuffers[ch];
        for (size_t i = 0; i < m_blockSize; ++i) {
            dst[i] = static_cast<double>(src[i]) * win[i];
        }
        
        // Check if we have enough frames for a batch
        size_t currentFrames = m_inputBuffer_ch[ch].size() / m_blockSize;
        if (currentFrames >= m_batchSize) {
            processBatch(ch, currentFrames);
            m_inputBuffer_ch[ch].clear();
        }
    }
    
    return fs;
}

void BIaccPlugin::processBatch(size_t channel, size_t numFrames)
{
    if (numFrames == 0) return;

    size_t blockSize = m_blockSize;
    double* accSpectrum = m_accumulatedSpectrum_ch[channel].data();
    float& globalMax = m_globalMax_ch[channel];
    
    for (size_t frame = 0; frame < numFrames; ++frame) {
        // Get pointer to current frame in input buffer
        double* frameData = &m_inputBuffer_ch[channel][frame * blockSize];
        
        // Perform FFT
        m_fft->forward(frameData, m_fftOut.data());
        const double* fftOut = m_fftOut.data();
        
        // Compute magnitudes and accumulate
        // Bin 0 (DC)
        double real0 = fftOut[0];
        double mag0 = std::abs(real0);
        if (static_cast<float>(mag0) > globalMax) globalMax = static_cast<float>(mag0);
        accSpectrum[0] += mag0 * mag0; // Accumulate Power

        // Bins 1 to N/2 - 1
        for (size_t i = 1; i < m_numBins; ++i) {
            double real = fftOut[2 * i];
            double imag = fftOut[2 * i + 1];
            double magSq = real * real + imag * imag;
            double magnitude = std::sqrt(magSq);
            
            if (static_cast<float>(magnitude) > globalMax) {
                globalMax = static_cast<float>(magnitude);
            }
            
            accSpectrum[i] += magSq;
        }
        
        m_frameCount_ch[channel]++;
    }
}

Vamp::Plugin::FeatureSet
BIaccPlugin::getRemainingFeatures()
{
    FeatureSet fs;
    
    // Process any remaining frames in the buffer for each channel
    for (size_t ch = 0; ch < m_channels; ++ch) {
        if (!m_inputBuffer_ch[ch].empty()) {
            size_t remainingFrames = m_inputBuffer_ch[ch].size() / m_blockSize;
            if (remainingFrames > 0) {
                processBatch(ch, remainingFrames);
            }
            m_inputBuffer_ch[ch].clear();
        }
    }
    
    // Calculate BI for each channel
    for (size_t ch = 0; ch < m_channels; ++ch) {
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
            if (val < 1e-12) val = 1e-12; // Floor to -120dB (Power) -> -120dB
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
