// ADIaccPlugin.cpp - Acoustic Diversity Index (Accumulated)
// Matches soundecology::acoustic_diversity behavior

#include "ADIaccPlugin.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <map>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

ADIaccPlugin::ADIaccPlugin(float inputSampleRate) :
    ACIBasePlugin(inputSampleRate),
    m_dbThreshold(-50.0f),
    m_freqStep(1000.0f)
{
    m_maxFreq = 10000.0f; // Default max freq for ADI
}

ADIaccPlugin::~ADIaccPlugin()
{
}

string
ADIaccPlugin::getIdentifier() const
{
    return "adi-acc";
}

string
ADIaccPlugin::getName() const
{
    return "Acoustic Diversity Index (Accumulated)";
}

string
ADIaccPlugin::getDescription() const
{
    return "Calculates the Acoustic Diversity Index (ADI) of a signal.";
}

string
ADIaccPlugin::getMaker() const
{
    return "ReVAMP";
}

int
ADIaccPlugin::getPluginVersion() const
{
    return 1;
}

string
ADIaccPlugin::getCopyright() const
{
    return "MIT License";
}

Vamp::Plugin::ParameterList
ADIaccPlugin::getParameterDescriptors() const
{
    ParameterList list;

    ParameterDescriptor d;
    d.identifier = "maxfreq";
    d.name = "Maximum Frequency";
    d.description = "Maximum frequency (kHz) to include in calculation";
    d.unit = "kHz";
    d.minValue = 0;
    d.maxValue = m_inputSampleRate / 2000.0f;
    d.defaultValue = 10.0f;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "db_threshold";
    d.name = "dB Threshold";
    d.description = "Threshold in dBFS";
    d.unit = "dB";
    d.minValue = -120.0f;
    d.maxValue = 0.0f;
    d.defaultValue = -50.0f;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "freq_step";
    d.name = "Frequency Step";
    d.description = "Size of frequency bands (Hz)";
    d.unit = "Hz";
    d.minValue = 100.0f;
    d.maxValue = 10000.0f;
    d.defaultValue = 1000.0f;
    d.isQuantized = false;
    list.push_back(d);

    return list;
}

float
ADIaccPlugin::getParameter(string identifier) const
{
    if (identifier == "maxfreq") return m_maxFreq; // Stored in base class but used here
    if (identifier == "db_threshold") return m_dbThreshold;
    if (identifier == "freq_step") return m_freqStep;
    return 0;
}

void
ADIaccPlugin::setParameter(string identifier, float value)
{
    if (identifier == "maxfreq") {
        m_maxFreq = value;
    } else if (identifier == "db_threshold") {
        m_dbThreshold = value;
    } else if (identifier == "freq_step") {
        m_freqStep = value;
    }
}

Vamp::Plugin::ProgramList
ADIaccPlugin::getPrograms() const
{
    ProgramList list;
    return list;
}

string
ADIaccPlugin::getCurrentProgram() const
{
    return "";
}

void
ADIaccPlugin::selectProgram(string name)
{
}

Vamp::Plugin::OutputList
ADIaccPlugin::getOutputDescriptors() const
{
    OutputList list;

    OutputDescriptor d;
    d.identifier = "adi";
    d.name = "ADI";
    d.description = "Acoustic Diversity Index";
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
ADIaccPlugin::getPreferredBlockSize() const
{
    return 4096;
}

size_t
ADIaccPlugin::getPreferredStepSize() const
{
    return 4096;
}

bool
ADIaccPlugin::initialise(size_t channels, size_t stepSize, size_t blockSize)
{
    if (!ACIBasePlugin::initialise(channels, stepSize, blockSize)) return false;

    // Reserve spectral data
    // ADI needs the whole file, so we reserve a reasonable amount
    m_spectralData.reserve(10000 * m_numBins);

    return true;
}

Vamp::Plugin::FeatureSet
ADIaccPlugin::getRemainingFeatures()
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
    
    if (m_spectralData.empty()) return fs;

    // 1. Find Global Max
    // soundecology normalizes to 0 dB.
    // This means finding the max magnitude and scaling so max = 1.0 (0 dB).
    
    size_t numFrames = m_spectralData.size() / m_numBins;
    size_t totalElements = numFrames * m_numBins;
    
    float globalMax = 0.0f;
    for (size_t i = 0; i < totalElements; ++i) {
        if (m_spectralData[i] > globalMax) {
            globalMax = m_spectralData[i];
        }
    }
    
    // 2. Optimization: Pre-calculate threshold and bin mapping
    // Instead of converting every bin to dB, we convert the threshold to linear magnitude.
    // dB = 20 * log10(mag / globalMax) > threshold
    // mag > globalMax * 10^(threshold / 20)
    
    float linearThreshold = 0.0f;
    if (globalMax > 0) {
        linearThreshold = globalMax * std::pow(10.0f, m_dbThreshold / 20.0f);
    } else {
        linearThreshold = 1e9f; // Unreachable if max is 0
    }
    
    // Pre-calculate bin to band mapping
    float maxFreqHz = m_maxFreq * 1000.0f; // m_maxFreq is in kHz
    if (maxFreqHz <= 0) maxFreqHz = m_inputSampleRate / 2.0f;
    
    int numBands = static_cast<int>(std::ceil(maxFreqHz / m_freqStep));
    std::vector<int> bandCounts(numBands, 0);
    int totalCount = 0;
    
    float binResolution = (m_inputSampleRate / 2.0f) / (m_blockSize / 2);
    
    std::vector<int> binToBand(m_numBins, -1);
    for (size_t b = 0; b < m_numBins; ++b) {
        float freq = (b + 1) * binResolution;
        if (freq <= maxFreqHz) {
            int bandIndex = static_cast<int>(std::floor(freq / m_freqStep));
            if (bandIndex >= 0 && bandIndex < numBands) {
                binToBand[b] = bandIndex;
            }
        }
    }
    
    // 3. Iterate and Count
    // Optimized loop: No log10, no division, no branching logic for bands inside loop
    for (size_t t = 0; t < numFrames; ++t) {
        size_t frameOffset = t * m_numBins;
        for (size_t b = 0; b < m_numBins; ++b) {
            if (m_spectralData[frameOffset + b] > linearThreshold) {
                int band = binToBand[b];
                if (band != -1) {
                    bandCounts[band]++;
                    totalCount++;
                }
            }
        }
    }
    
    // 4. Calculate Shannon Index
    // H = - sum(p * ln(p))
    double shannon = 0.0;
    if (totalCount > 0) {
        for (int count : bandCounts) {
            if (count > 0) {
                double p = static_cast<double>(count) / totalCount;
                shannon -= p * std::log(p);
            }
        }
    }
    
    Feature f;
    f.hasTimestamp = true;
    f.timestamp = Vamp::RealTime::frame2RealTime(m_frameCount * m_stepSize, m_inputSampleRate);
    f.values.push_back(static_cast<float>(shannon));
    
    fs[0].push_back(f);
    
    return fs;
}
