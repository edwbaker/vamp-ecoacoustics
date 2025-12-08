// BIaccPlugin.cpp - Bioacoustic Index (Accumulated)
// Matches soundecology::bioacoustic_index behavior

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
    EcoacousticSpectralPlugin(inputSampleRate),
    m_minFreq(2000.0f),
    m_maxFreq(8000.0f)
{
}

BIaccPlugin::~BIaccPlugin()
{
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
    return "ReVAMP";
}

int
BIaccPlugin::getPluginVersion() const
{
    return 1;
}

string
BIaccPlugin::getCopyright() const
{
    return "MIT License";
}

Vamp::Plugin::ParameterList
BIaccPlugin::getParameterDescriptors() const
{
    ParameterList list;

    ParameterDescriptor d;
    d.identifier = "min_freq";
    d.name = "Minimum Frequency";
    d.description = "Minimum frequency (Hz)";
    d.unit = "Hz";
    d.minValue = 0;
    d.maxValue = m_inputSampleRate / 2.0f;
    d.defaultValue = 2000.0f;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "max_freq";
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
    if (identifier == "min_freq") return m_minFreq;
    if (identifier == "max_freq") return m_maxFreq;
    return 0;
}

void
BIaccPlugin::setParameter(string identifier, float value)
{
    if (identifier == "min_freq") {
        m_minFreq = value;
    } else if (identifier == "max_freq") {
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
}

Vamp::Plugin::OutputList
BIaccPlugin::getOutputDescriptors() const
{
    OutputList list;

    OutputDescriptor d;
    d.identifier = "bi";
    d.name = "Bioacoustic Index";
    d.description = "Bioacoustic Index";
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
BIaccPlugin::getPreferredBlockSize() const
{
    return 512; // soundecology default
}

size_t
BIaccPlugin::getPreferredStepSize() const
{
    return 512; // No overlap default
}

bool
BIaccPlugin::initialise(size_t channels, size_t stepSize, size_t blockSize)
{
    if (!EcoacousticSpectralPlugin::initialise(channels, stepSize, blockSize)) return false;

    m_accumulatedSpectrum.assign(m_numBins, 0.0);

    return true;
}

void BIaccPlugin::processBatch(size_t numFrames)
{
    if (numFrames == 0) return;

    size_t blockSize = m_blockSize;
    
    for (size_t frame = 0; frame < numFrames; ++frame) {
        // Get pointer to current frame in input buffer
        double* frameData = &m_inputBuffer[frame * blockSize];
        
        // Perform FFT
        m_fft->forward(frameData, m_fftOut.data());
        
        // Compute magnitudes and accumulate
        for (size_t i = 1; i <= blockSize / 2; ++i) {
            double real = m_fftOut[2 * i];
            double imag = m_fftOut[2 * i + 1];
            double magnitude = std::sqrt(real * real + imag * imag);
            
            if (magnitude > m_globalMax) {
                m_globalMax = magnitude;
            }
            
            m_accumulatedSpectrum[i-1] += magnitude;
        }
        
        m_frameCount++;
    }
}

Vamp::Plugin::FeatureSet
BIaccPlugin::getRemainingFeatures()
{
    FeatureSet fs;
    
    // Process any remaining frames in the buffer
    if (!m_inputBuffer.empty()) {
        size_t remainingFrames = m_inputBuffer.size() / m_blockSize;
        if (remainingFrames > 0) {
            processBatch(remainingFrames);
        }
        m_inputBuffer.clear();
    }
    
    if (m_frameCount == 0) return fs;

    // 1. Calculate Mean Spectrum
    std::vector<double> meanSpectrum(m_numBins);
    
    for (size_t b = 0; b < m_numBins; ++b) {
        meanSpectrum[b] = m_accumulatedSpectrum[b] / m_frameCount;
    }
    
    // 2. Convert to dB Normalized (Max 0)
    // soundecology logic:
    // spec_left <- spectro(..., dB="max0")$amp  -> Returns dB normalized to max=0
    // specA_left <- apply(spec_left, 1, meandB) -> Mean of dB (effectively 20*log10(Mean(Amp)/MaxAmp))
    
    float maxAmp = m_globalMax;
    if (maxAmp <= 1e-9f) maxAmp = 1.0f; 

    std::vector<double> meanSpectrumDB(m_numBins);
    for (size_t b = 0; b < m_numBins; ++b) {
        double val = meanSpectrum[b] / maxAmp;
        if (val < 1e-12) val = 1e-12; // Floor to -240dB
        meanSpectrumDB[b] = 20.0 * std::log10(val);
    }

    // 3. Extract Segment and Normalize
    // rows_width = length(specA_left)/nyquist_freq = (N/2) / (fs/2) = N/fs = 1/df
    float df = (m_inputSampleRate / 2.0f) / m_numBins;
    float rows_width = 1.0f / df;
    
    // R uses 1-based indexing, we use 0-based.
    // min_row = min_freq * rows_width
    size_t minBin = static_cast<size_t>(std::floor(m_minFreq * rows_width));
    size_t maxBin = static_cast<size_t>(std::ceil(m_maxFreq * rows_width));
    
    if (minBin >= m_numBins) minBin = m_numBins - 1;
    if (maxBin > m_numBins) maxBin = m_numBins;
    
    if (minBin >= maxBin) {
        Feature f;
        f.hasTimestamp = true;
        f.timestamp = Vamp::RealTime::frame2RealTime(0, m_inputSampleRate);
        f.values.push_back(0.0f);
        fs[0].push_back(f);
        return fs;
    }

    // Find min in segment
    double minVal = 1e9;
    for (size_t b = minBin; b < maxBin; ++b) {
        if (meanSpectrumDB[b] < minVal) minVal = meanSpectrumDB[b];
    }
    
    // 4. Calculate Area
    // left_area <- sum((specA_left_segment - min) * rows_width)
    double bi = 0.0;
    for (size_t b = minBin; b < maxBin; ++b) {
        bi += (meanSpectrumDB[b] - minVal) * rows_width;
    }
    
    Feature f;
    f.hasTimestamp = true;
    f.timestamp = Vamp::RealTime::frame2RealTime(0, m_inputSampleRate);
    f.values.push_back(static_cast<float>(bi));
    
    fs[0].push_back(f);
    
    return fs;
}
