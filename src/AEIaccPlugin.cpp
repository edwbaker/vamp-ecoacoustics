// AEIaccPlugin.cpp - Acoustic Evenness Index (Accumulated)
// Matches soundecology::acoustic_eveness behavior

#include "AEIaccPlugin.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

AEIaccPlugin::AEIaccPlugin(float inputSampleRate) :
    ACIBasePlugin(inputSampleRate),
    m_dbThreshold(-50.0f),
    m_freqStep(1000.0f)
{
    m_maxFreq = 10000.0f; // Default max freq for AEI in soundecology
}

AEIaccPlugin::~AEIaccPlugin()
{
}

string
AEIaccPlugin::getIdentifier() const
{
    return "aei-acc";
}

string
AEIaccPlugin::getName() const
{
    return "Acoustic Evenness Index (soundecology)";
}

string
AEIaccPlugin::getDescription() const
{
    return "Calculates the Acoustic Evenness Index (AEI) of a signal, matching soundecology implementation.";
}

string
AEIaccPlugin::getMaker() const
{
    return "ReVAMP";
}

int
AEIaccPlugin::getPluginVersion() const
{
    return 1;
}

string
AEIaccPlugin::getCopyright() const
{
    return "MIT License";
}

Vamp::Plugin::ParameterList
AEIaccPlugin::getParameterDescriptors() const
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
AEIaccPlugin::getParameter(string identifier) const
{
    if (identifier == "maxfreq") return m_maxFreq / 1000.0f;
    if (identifier == "db_threshold") return m_dbThreshold;
    if (identifier == "freq_step") return m_freqStep;
    return 0;
}

void
AEIaccPlugin::setParameter(string identifier, float value)
{
    if (identifier == "maxfreq") {
        m_maxFreq = value * 1000.0f;
    } else if (identifier == "db_threshold") {
        m_dbThreshold = value;
    } else if (identifier == "freq_step") {
        m_freqStep = value;
    }
}

Vamp::Plugin::ProgramList
AEIaccPlugin::getPrograms() const
{
    ProgramList list;
    return list;
}

string
AEIaccPlugin::getCurrentProgram() const
{
    return "";
}

void
AEIaccPlugin::selectProgram(string name)
{
}

Vamp::Plugin::OutputList
AEIaccPlugin::getOutputDescriptors() const
{
    OutputList list;

    OutputDescriptor d;
    d.identifier = "aei";
    d.name = "AEI";
    d.description = "Acoustic Evenness Index";
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
AEIaccPlugin::getPreferredBlockSize() const
{
    return 4096; // soundecology default window length is usually sample_rate/10, but we use fixed block size
}

size_t
AEIaccPlugin::getPreferredStepSize() const
{
    return 4096; // No overlap by default in soundecology?
}

bool
AEIaccPlugin::initialise(size_t channels, size_t stepSize, size_t blockSize)
{
    if (!ACIBasePlugin::initialise(channels, stepSize, blockSize)) return false;

    // Reserve spectral data
    m_spectralData.reserve(10000 * m_numBins);

    return true;
}

Vamp::Plugin::FeatureSet
AEIaccPlugin::getRemainingFeatures()
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
    // Optimized: m_globalMax is tracked in ACIBasePlugin
    size_t numFrames = m_spectralData.size() / m_numBins;
    float globalMax = m_globalMax;
    
    // 2. Calculate Threshold
    float linearThreshold = 0.0f;
    if (globalMax > 0) {
        linearThreshold = globalMax * std::pow(10.0f, m_dbThreshold / 20.0f);
    } else {
        linearThreshold = 1e9f; // Unreachable
    }
    
    // 3. Pre-calculate bin to band mapping
    float maxFreqHz = m_maxFreq; // m_maxFreq is already in Hz (stored in ACIBasePlugin)
    // Wait, setParameter sets it in Hz? 
    // In ADIaccPlugin setParameter: m_maxFreq = value * 1000.0f;
    // So m_maxFreq is in Hz.
    
    if (maxFreqHz <= 0) maxFreqHz = m_inputSampleRate / 2.0f;
    
    int numBands = static_cast<int>(std::ceil(maxFreqHz / m_freqStep));
    // soundecology usually has 10 bands of 1000Hz up to 10000Hz.
    
    std::vector<double> bandCounts(numBands, 0.0); // Using double for Gini calc
    
    float binResolution = (m_inputSampleRate / 2.0f) / (m_blockSize / 2);
    
    struct BandRange {
        int startBin;
        int endBin; // Exclusive
    };
    std::vector<BandRange> bandRanges(numBands, {-1, -1});

    for (size_t b = 0; b < m_numBins; ++b) {
        float freq = (b + 1) * binResolution;
        if (freq <= maxFreqHz) {
            int bandIndex = static_cast<int>(std::floor(freq / m_freqStep));
            if (bandIndex >= 0 && bandIndex < numBands) {
                if (bandRanges[bandIndex].startBin == -1) {
                    bandRanges[bandIndex].startBin = b;
                }
                bandRanges[bandIndex].endBin = b + 1;
            }
        }
    }
    
    // 4. Iterate and Count
    for (size_t t = 0; t < numFrames; ++t) {
        size_t frameOffset = t * m_numBins;
        
        for (int band = 0; band < numBands; ++band) {
            int start = bandRanges[band].startBin;
            int end = bandRanges[band].endBin;
            
            if (start != -1 && end != -1) {
                for (int b = start; b < end; ++b) {
                    if (m_spectralData[frameOffset + b] > linearThreshold) {
                        bandCounts[band] += 1.0;
                    }
                }
            }
        }
    }
    
    // 5. Calculate Gini Index
    // G = (2 * sum(i * xi) / (n * sum(xi))) - (n + 1) / n
    // where xi are sorted in ascending order, i is 1-based index
    
    std::sort(bandCounts.begin(), bandCounts.end());
    
    double sumCounts = 0.0;
    double weightedSum = 0.0;
    double n = static_cast<double>(numBands);
    
    for (int i = 0; i < numBands; ++i) {
        sumCounts += bandCounts[i];
        weightedSum += (i + 1) * bandCounts[i];
    }
    
    double aei = 0.0;
    if (sumCounts > 0) {
        aei = (2.0 * weightedSum) / (n * sumCounts) - (n + 1.0) / n;
    } else {
        // If all counts are 0, Gini is 0 (perfect equality of nothing?)
        // Or undefined. soundecology returns 0 if no signal?
        aei = 0.0;
    }
    
    Feature f;
    f.hasTimestamp = true;
    f.timestamp = Vamp::RealTime::frame2RealTime(m_frameCount * m_stepSize, m_inputSampleRate);
    f.values.push_back(static_cast<float>(aei));
    
    fs[0].push_back(f);
    
    return fs;
}
