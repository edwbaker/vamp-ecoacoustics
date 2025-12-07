// AEIPlugin.cpp - Acoustic Evenness Index
// Matches scikit-maad acoustic_eveness_index behavior

#include "AEIPlugin.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

AEIPlugin::AEIPlugin(float inputSampleRate) :
    ACIBasePlugin(inputSampleRate),
    m_minFreq(0.0f),
    m_maxFreq(20000.0f),
    m_binStep(1000.0f),
    m_dbThreshold(-50.0f)
{
}

AEIPlugin::~AEIPlugin()
{
}

string
AEIPlugin::getIdentifier() const
{
    return "aei";
}

string
AEIPlugin::getName() const
{
    return "Acoustic Evenness Index (scikit-maad)";
}

string
AEIPlugin::getDescription() const
{
    return "Calculates the Acoustic Evenness Index (AEI) of a signal, matching scikit-maad implementation.";
}

string
AEIPlugin::getMaker() const
{
    return "ReVAMP";
}

int
AEIPlugin::getPluginVersion() const
{
    return 1;
}

string
AEIPlugin::getCopyright() const
{
    return "MIT License";
}

Vamp::Plugin::ParameterList
AEIPlugin::getParameterDescriptors() const
{
    ParameterList list;

    ParameterDescriptor d;
    d.identifier = "min_freq";
    d.name = "Minimum Frequency";
    d.description = "Minimum frequency (Hz) to include in calculation";
    d.unit = "Hz";
    d.minValue = 0;
    d.maxValue = m_inputSampleRate / 2.0f;
    d.defaultValue = 0.0f;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "max_freq";
    d.name = "Maximum Frequency";
    d.description = "Maximum frequency (Hz) to include in calculation";
    d.unit = "Hz";
    d.minValue = 0;
    d.maxValue = m_inputSampleRate / 2.0f;
    d.defaultValue = 20000.0f;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "bin_step";
    d.name = "Frequency Step";
    d.description = "Size of frequency bands (Hz)";
    d.unit = "Hz";
    d.minValue = 100.0f;
    d.maxValue = 10000.0f;
    d.defaultValue = 1000.0f;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "db_threshold";
    d.name = "dB Threshold";
    d.description = "Threshold in dBFS (relative to global max)";
    d.unit = "dB";
    d.minValue = -120.0f;
    d.maxValue = 0.0f;
    d.defaultValue = -50.0f;
    d.isQuantized = false;
    list.push_back(d);

    return list;
}

float
AEIPlugin::getParameter(string identifier) const
{
    if (identifier == "min_freq") return m_minFreq;
    if (identifier == "max_freq") return m_maxFreq;
    if (identifier == "bin_step") return m_binStep;
    if (identifier == "db_threshold") return m_dbThreshold;
    return 0;
}

void
AEIPlugin::setParameter(string identifier, float value)
{
    if (identifier == "min_freq") {
        m_minFreq = value;
    } else if (identifier == "max_freq") {
        m_maxFreq = value;
    } else if (identifier == "bin_step") {
        m_binStep = value;
    } else if (identifier == "db_threshold") {
        m_dbThreshold = value;
    }
}

Vamp::Plugin::ProgramList
AEIPlugin::getPrograms() const
{
    ProgramList list;
    return list;
}

string
AEIPlugin::getCurrentProgram() const
{
    return "";
}

void
AEIPlugin::selectProgram(string name)
{
}

Vamp::Plugin::OutputList
AEIPlugin::getOutputDescriptors() const
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
AEIPlugin::getPreferredBlockSize() const
{
    return 4096;
}

size_t
AEIPlugin::getPreferredStepSize() const
{
    return 4096;
}

bool
AEIPlugin::initialise(size_t channels, size_t stepSize, size_t blockSize)
{
    if (!ACIBasePlugin::initialise(channels, stepSize, blockSize)) return false;

    // Reserve spectral data
    m_spectralData.reserve(10000 * m_numBins);

    return true;
}

Vamp::Plugin::FeatureSet
AEIPlugin::getRemainingFeatures()
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
    
    // 3. Calculate Scores per Band
    int numBands = static_cast<int>(std::floor((m_maxFreq - m_minFreq) / m_binStep));
    if (numBands < 1) numBands = 1;

    std::vector<double> scores(numBands, 0.0);
    std::vector<int> binCounts(numBands, 0); // Total bins in each band (across all time)
    std::vector<int> hitCounts(numBands, 0); // Bins > threshold
    
    float binResolution = (m_inputSampleRate / 2.0f) / (m_blockSize / 2);
    
    // Map bands to bin ranges
    // Instead of per-bin lookup, we precompute which bins belong to which band.
    // Since bands are contiguous frequency ranges, we can store start/end bin indices.
    struct BandRange {
        int startBin;
        int endBin; // Exclusive
    };
    std::vector<BandRange> bandRanges(numBands, {-1, -1});

    for (size_t b = 0; b < m_numBins; ++b) {
        float freq = (b + 1) * binResolution;
        
        if (freq >= m_minFreq && freq < m_maxFreq) {
            int bandIndex = static_cast<int>(std::floor((freq - m_minFreq) / m_binStep));
            if (bandIndex >= 0 && bandIndex < numBands) {
                if (bandRanges[bandIndex].startBin == -1) {
                    bandRanges[bandIndex].startBin = b;
                }
                bandRanges[bandIndex].endBin = b + 1;
            }
        }
    }
    
    // Iterate
    // Optimized loop: Iterate frames, then bands, then bins in range.
    // This avoids the indirect lookup binToBand[b] and the check if band != -1.
    for (size_t t = 0; t < numFrames; ++t) {
        size_t frameOffset = t * m_numBins;
        
        for (int band = 0; band < numBands; ++band) {
            int start = bandRanges[band].startBin;
            int end = bandRanges[band].endBin;
            
            if (start != -1 && end != -1) {
                int count = end - start;
                binCounts[band] += count;
                
                for (int b = start; b < end; ++b) {
                    if (m_spectralData[frameOffset + b] > linearThreshold) {
                        hitCounts[band]++;
                    }
                }
            }
        }
    }
    
    // Calculate scores (proportions)
    for (int i = 0; i < numBands; ++i) {
        if (binCounts[i] > 0) {
            scores[i] = static_cast<double>(hitCounts[i]) / binCounts[i];
        } else {
            scores[i] = 0.0;
        }
    }
    
    // 4. Calculate Gini Index
    // G = (2 * sum(i * xi) / (n * sum(xi))) - (n + 1) / n
    // where xi are sorted in ascending order, i is 1-based index
    
    std::sort(scores.begin(), scores.end());
    
    double sumScores = 0.0;
    double weightedSum = 0.0;
    double n = static_cast<double>(numBands);
    
    for (int i = 0; i < numBands; ++i) {
        sumScores += scores[i];
        weightedSum += (i + 1) * scores[i];
    }
    
    double aei = 0.0;
    if (sumScores > 0) {
        aei = (2.0 * weightedSum) / (n * sumScores) - (n + 1.0) / n;
    } else {
        // If sum is 0, all scores are 0. Perfect equality?
        // Gini of [0, 0, 0] is technically undefined (0/0) or 0.
        // scikit-maad likely returns 0 or NaN.
        // Let's assume 0 (perfect evenness of silence).
        aei = 0.0;
    }
    
    Feature f;
    f.hasTimestamp = true;
    f.timestamp = Vamp::RealTime::frame2RealTime(m_frameCount * m_stepSize, m_inputSampleRate);
    f.values.push_back(static_cast<float>(aei));
    
    fs[0].push_back(f);
    
    return fs;
}
