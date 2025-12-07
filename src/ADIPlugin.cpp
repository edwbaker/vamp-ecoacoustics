// ADIPlugin.cpp - Acoustic Diversity Index
// Matches scikit-maad acoustic_diversity_index behavior

#include "ADIPlugin.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

ADIPlugin::ADIPlugin(float inputSampleRate) :
    ACIBasePlugin(inputSampleRate),
    m_minFreq(0.0f),
    m_maxFreq(20000.0f),
    m_binStep(1000.0f),
    m_dbThreshold(-50.0f),
    m_indexType(0) // Shannon
{
}

ADIPlugin::~ADIPlugin()
{
}

string
ADIPlugin::getIdentifier() const
{
    return "adi";
}

string
ADIPlugin::getName() const
{
    return "Acoustic Diversity Index (scikit-maad)";
}

string
ADIPlugin::getDescription() const
{
    return "Calculates the Acoustic Diversity Index (ADI) of a signal, matching scikit-maad implementation.";
}

string
ADIPlugin::getMaker() const
{
    return "ReVAMP";
}

int
ADIPlugin::getPluginVersion() const
{
    return 1;
}

string
ADIPlugin::getCopyright() const
{
    return "MIT License";
}

Vamp::Plugin::ParameterList
ADIPlugin::getParameterDescriptors() const
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

    d.identifier = "index_type";
    d.name = "Index Type";
    d.description = "Diversity Index Type (0: Shannon, 1: Simpson, 2: Inverse Simpson)";
    d.unit = "";
    d.minValue = 0;
    d.maxValue = 2;
    d.defaultValue = 0;
    d.isQuantized = true;
    d.quantizeStep = 1.0f;
    d.valueNames.push_back("Shannon");
    d.valueNames.push_back("Simpson");
    d.valueNames.push_back("Inverse Simpson");
    list.push_back(d);

    return list;
}

float
ADIPlugin::getParameter(string identifier) const
{
    if (identifier == "min_freq") return m_minFreq;
    if (identifier == "max_freq") return m_maxFreq;
    if (identifier == "bin_step") return m_binStep;
    if (identifier == "db_threshold") return m_dbThreshold;
    if (identifier == "index_type") return (float)m_indexType;
    return 0;
}

void
ADIPlugin::setParameter(string identifier, float value)
{
    if (identifier == "min_freq") {
        m_minFreq = value;
    } else if (identifier == "max_freq") {
        m_maxFreq = value;
    } else if (identifier == "bin_step") {
        m_binStep = value;
    } else if (identifier == "db_threshold") {
        m_dbThreshold = value;
    } else if (identifier == "index_type") {
        m_indexType = (int)value;
    }
}

Vamp::Plugin::ProgramList
ADIPlugin::getPrograms() const
{
    ProgramList list;
    return list;
}

string
ADIPlugin::getCurrentProgram() const
{
    return "";
}

void
ADIPlugin::selectProgram(string name)
{
}

Vamp::Plugin::OutputList
ADIPlugin::getOutputDescriptors() const
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
ADIPlugin::getPreferredBlockSize() const
{
    return 4096;
}

size_t
ADIPlugin::getPreferredStepSize() const
{
    return 4096;
}

bool
ADIPlugin::initialise(size_t channels, size_t stepSize, size_t blockSize)
{
    if (!ACIBasePlugin::initialise(channels, stepSize, blockSize)) return false;

    // Reserve spectral data
    m_spectralData.reserve(10000 * m_numBins);

    return true;
}

Vamp::Plugin::FeatureSet
ADIPlugin::getRemainingFeatures()
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
    size_t numFrames = m_spectralData.size() / m_numBins;
    size_t totalElements = numFrames * m_numBins;
    
    float globalMax = 0.0f;
    for (size_t i = 0; i < totalElements; ++i) {
        if (m_spectralData[i] > globalMax) {
            globalMax = m_spectralData[i];
        }
    }
    
    // 2. Calculate Threshold
    // scikit-maad: Sxx_dB = amplitude2dB(Sxx/max(Sxx))
    // threshold is in dB.
    // val_dB = 20 * log10(val / globalMax)
    // val > globalMax * 10^(threshold/20)
    
    float linearThreshold = 0.0f;
    if (globalMax > 0) {
        linearThreshold = globalMax * std::pow(10.0f, m_dbThreshold / 20.0f);
    } else {
        linearThreshold = 1e9f; // Unreachable
    }
    
    // 3. Calculate Scores per Band
    // N = floor((fmax-fmin)/bin_step)
    int numBands = static_cast<int>(std::floor((m_maxFreq - m_minFreq) / m_binStep));
    if (numBands < 1) numBands = 1;

    std::vector<double> scores(numBands, 0.0);
    std::vector<int> binCounts(numBands, 0); // Total bins in each band (across all time)
    std::vector<int> hitCounts(numBands, 0); // Bins > threshold
    
    float binResolution = (m_inputSampleRate / 2.0f) / (m_blockSize / 2);
    
    // Map bins to bands
    std::vector<int> binToBand(m_numBins, -1);
    for (size_t b = 0; b < m_numBins; ++b) {
        float freq = b * binResolution; // Center freq? Or start? scikit-maad uses fn vector.
        // Usually bin 0 is DC (0 Hz).
        // scikit-maad: f0 = fmin + step*ii; f1 = f0 + step.
        // index_bw selects bins where freq >= f0 and freq <= f1 (inclusive usually in maad/numpy logic if using boolean mask)
        // But we should be careful about overlap.
        // Let's use [f0, f1).
        
        if (freq >= m_minFreq && freq < m_maxFreq) {
            int bandIndex = static_cast<int>(std::floor((freq - m_minFreq) / m_binStep));
            if (bandIndex >= 0 && bandIndex < numBands) {
                binToBand[b] = bandIndex;
            }
        }
    }
    
    // Iterate
    for (size_t t = 0; t < numFrames; ++t) {
        size_t frameOffset = t * m_numBins;
        for (size_t b = 0; b < m_numBins; ++b) {
            int band = binToBand[b];
            if (band != -1) {
                binCounts[band]++;
                if (m_spectralData[frameOffset + b] > linearThreshold) {
                    hitCounts[band]++;
                }
            }
        }
    }
    
    // Calculate scores (proportions)
    double sumScores = 0.0;
    for (int i = 0; i < numBands; ++i) {
        if (binCounts[i] > 0) {
            scores[i] = static_cast<double>(hitCounts[i]) / binCounts[i];
        } else {
            scores[i] = 0.0;
        }
        sumScores += scores[i];
    }
    
    // 4. Calculate Index
    double adi = 0.0;
    
    if (sumScores > 0) {
        // Normalize scores to probability distribution
        std::vector<double> p(numBands);
        for (int i = 0; i < numBands; ++i) {
            p[i] = scores[i] / sumScores;
        }
        
        if (m_indexType == 0) { // Shannon
            for (double val : p) {
                if (val > 0) {
                    adi -= val * std::log(val);
                }
            }
        } else if (m_indexType == 1) { // Simpson
            double sumSq = 0.0;
            for (double val : p) {
                sumSq += val * val;
            }
            adi = 1.0 - sumSq;
        } else if (m_indexType == 2) { // Inverse Simpson
            double sumSq = 0.0;
            for (double val : p) {
                sumSq += val * val;
            }
            if (sumSq > 0) {
                adi = 1.0 / sumSq;
            }
        }
    }
    
    Feature f;
    f.hasTimestamp = true;
    f.timestamp = Vamp::RealTime::frame2RealTime(m_frameCount * m_stepSize, m_inputSampleRate);
    f.values.push_back(static_cast<float>(adi));
    
    fs[0].push_back(f);
    
    return fs;
}
