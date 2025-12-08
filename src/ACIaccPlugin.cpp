// ACIaccPlugin.cpp - Acoustic Complexity Index (Accumulated)
// Matches soundecology::acoustic_complexity behavior

#include "ACIaccPlugin.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

ACIaccPlugin::ACIaccPlugin(float inputSampleRate) :
    EcoacousticSpectralPlugin(inputSampleRate),
    m_clusterSize(5.0f), // Default 5 seconds
    m_runningTotalACI(0.0f),
    m_framesPerCluster(0)
{
}

ACIaccPlugin::~ACIaccPlugin()
{
}

string
ACIaccPlugin::getIdentifier() const
{
    return "aci-acc";
}

string
ACIaccPlugin::getName() const
{
    return "Acoustic Complexity Index (Accumulated)";
}

string
ACIaccPlugin::getDescription() const
{
    return "Calculates ACI matching soundecology::acoustic_complexity (accumulated over clusters).";
}

string
ACIaccPlugin::getMaker() const
{
    return "ReVAMP";
}

int
ACIaccPlugin::getPluginVersion() const
{
    return 1;
}

string
ACIaccPlugin::getCopyright() const
{
    return "MIT License";
}

Vamp::Plugin::ParameterList
ACIaccPlugin::getParameterDescriptors() const
{
    ParameterList list;

    ParameterDescriptor d;
    d.identifier = "minfreq";
    d.name = "Minimum Frequency";
    d.description = "Minimum frequency (kHz) to include in calculation";
    d.unit = "kHz";
    d.minValue = 0;
    d.maxValue = m_inputSampleRate / 2000.0f;
    d.defaultValue = 0;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "maxfreq";
    d.name = "Maximum Frequency";
    d.description = "Maximum frequency (kHz) to include in calculation";
    d.unit = "kHz";
    d.minValue = 0;
    d.maxValue = m_inputSampleRate / 2000.0f;
    d.defaultValue = 0; // 0 means max
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "clustersize";
    d.name = "Cluster Size";
    d.description = "Size of the cluster in seconds (j)";
    d.unit = "s";
    d.minValue = 0.1;
    d.maxValue = 600;
    d.defaultValue = 5.0;
    d.isQuantized = false;
    list.push_back(d);

    return list;
}

float
ACIaccPlugin::getParameter(string identifier) const
{
    if (identifier == "minfreq") return m_minFreq;
    if (identifier == "maxfreq") return m_maxFreq;
    if (identifier == "clustersize") return m_clusterSize;
    return 0;
}

void
ACIaccPlugin::setParameter(string identifier, float value)
{
    if (identifier == "minfreq") {
        m_minFreq = value;
    } else if (identifier == "maxfreq") {
        m_maxFreq = value;
    } else if (identifier == "clustersize") {
        m_clusterSize = value;
        // Update frames per cluster if we are already initialized
        if (m_stepSize > 0 && m_inputSampleRate > 0) {
             float framesPerSecond = m_inputSampleRate / m_stepSize;
             m_framesPerCluster = static_cast<size_t>(std::floor(m_clusterSize * framesPerSecond));
             if (m_framesPerCluster < 1) m_framesPerCluster = 1;
        }
    }
}

Vamp::Plugin::ProgramList
ACIaccPlugin::getPrograms() const
{
    ProgramList list;
    return list;
}

string
ACIaccPlugin::getCurrentProgram() const
{
    return "";
}

void
ACIaccPlugin::selectProgram(string name)
{
}

Vamp::Plugin::OutputList
ACIaccPlugin::getOutputDescriptors() const
{
    OutputList list;

    OutputDescriptor d;
    d.identifier = "aci";
    d.name = "ACI";
    d.description = "Acoustic Complexity Index total value";
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

bool
ACIaccPlugin::initialise(size_t channels, size_t stepSize, size_t blockSize)
{
    if (!EcoacousticSpectralPlugin::initialise(channels, stepSize, blockSize)) return false;

    // Calculate frames per cluster
    float framesPerSecond = m_inputSampleRate / m_stepSize;
    m_framesPerCluster = static_cast<size_t>(std::floor(m_clusterSize * framesPerSecond));
    if (m_framesPerCluster < 1) m_framesPerCluster = 1;

    // Reserve spectral data for at least one cluster + batch
    m_spectralData.reserve((m_framesPerCluster + m_batchSize) * m_numBins);
    
    // Resize accumulators
    m_sumIntensity.resize(m_numBins);
    m_sumAbsDiff.resize(m_numBins);

    return true;
}

void
ACIaccPlugin::reset()
{
    EcoacousticSpectralPlugin::reset();
    m_runningTotalACI = 0.0f;
}

void ACIaccPlugin::processBatch(size_t numFrames)
{
    if (numFrames == 0) return;

    // Compute magnitudes using base class helper
    computeMagnitudes(numFrames);

    // Check if we have enough frames for a cluster
    // m_spectralData contains flattened frames.
    // Number of frames currently stored:
    size_t storedFrames = m_spectralData.size() / m_numBins;

    while (storedFrames >= m_framesPerCluster) {
        // Process one cluster
        size_t startFrame = 0;
        size_t endFrame = m_framesPerCluster;
        
        // Determine frequency range
        size_t minBinIndex = 0;
        size_t maxBinIndex = m_numBins - 1;
        
        if (m_minFreq > 0 || m_maxFreq > 0) {
            float binResolution = (m_inputSampleRate / 2.0f) / (m_blockSize / 2);
            if (m_minFreq > 0) {
                minBinIndex = static_cast<size_t>(std::round(m_minFreq * 1000.0f / binResolution));
            }
            if (m_maxFreq > 0) {
                maxBinIndex = static_cast<size_t>(std::round(m_maxFreq * 1000.0f / binResolution));
                if (maxBinIndex >= m_numBins) maxBinIndex = m_numBins - 1;
            }
        }

        size_t numActiveBins = maxBinIndex - minBinIndex + 1;
        
        // Reset accumulators
        std::fill(m_sumIntensity.begin(), m_sumIntensity.end(), 0.0f);
        std::fill(m_sumAbsDiff.begin(), m_sumAbsDiff.end(), 0.0f);
        
        // First frame of the cluster
        size_t firstFrameIdx = startFrame * m_numBins;
        for (size_t b = 0; b < numActiveBins; ++b) {
            m_sumIntensity[b] += m_spectralData[firstFrameIdx + minBinIndex + b];
        }
        
        // Subsequent frames
        for (size_t t = startFrame; t < endFrame - 1; ++t) {
            size_t currentFrameIdx = t * m_numBins;
            size_t nextFrameIdx = (t + 1) * m_numBins;
            
            for (size_t b = 0; b < numActiveBins; ++b) {
                float valCurrent = m_spectralData[currentFrameIdx + minBinIndex + b];
                float valNext = m_spectralData[nextFrameIdx + minBinIndex + b];
                
                m_sumIntensity[b] += valNext;
                m_sumAbsDiff[b] += std::abs(valNext - valCurrent);
            }
        }
        
        // Calculate ACI for this cluster
        float clusterAci = 0.0f;
        for (size_t b = 0; b < numActiveBins; ++b) {
            if (m_sumIntensity[b] > 0) {
                clusterAci += m_sumAbsDiff[b] / m_sumIntensity[b];
            }
        }
        
        m_runningTotalACI += clusterAci;

        // Remove processed frames
        // This is O(N) where N is remaining data size. 
        // Since we process as soon as we have a cluster, remaining data is small (< 1 cluster).
        // So this is efficient.
        size_t elementsToRemove = m_framesPerCluster * m_numBins;
        m_spectralData.erase(m_spectralData.begin(), m_spectralData.begin() + elementsToRemove);
        
        storedFrames -= m_framesPerCluster;
    }
}

Vamp::Plugin::FeatureSet
ACIaccPlugin::getRemainingFeatures()
{
    FeatureSet fs;
    
    // Process remaining frames
    if (!m_inputBuffer.empty()) {
        size_t remainingFrames = m_inputBuffer.size() / m_blockSize;
        if (remainingFrames > 0) {
            processBatch(remainingFrames);
        }
        m_inputBuffer.clear();
    }
    
    // Any remaining spectral data in m_spectralData is less than one cluster
    // soundecology ignores partial clusters at the end (floor division)
    // So we just return the accumulated total.

    Feature f;
    f.hasTimestamp = false;
    f.values.push_back(m_runningTotalACI);
    
    fs[0].push_back(f);
    
    return fs;
}
