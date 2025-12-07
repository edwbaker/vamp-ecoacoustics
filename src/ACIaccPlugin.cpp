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
    Plugin(inputSampleRate),
    m_minFreq(0),
    m_maxFreq(0),
    m_clusterSize(5.0f), // Default 5 seconds
    m_channels(0),
    m_blockSize(0),
    m_stepSize(0),
    m_numBins(0),
    m_fft(nullptr),
    m_batchSize(256),
    m_frameCount(0)
{
}

ACIaccPlugin::~ACIaccPlugin()
{
    if (m_fft) {
        delete m_fft;
        m_fft = nullptr;
    }
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

Vamp::Plugin::InputDomain
ACIaccPlugin::getInputDomain() const
{
    return TimeDomain;
}

size_t
ACIaccPlugin::getPreferredBlockSize() const
{
    return 512;
}

size_t
ACIaccPlugin::getPreferredStepSize() const
{
    return 512;
}

size_t
ACIaccPlugin::getMinChannelCount() const
{
    return 1;
}

size_t
ACIaccPlugin::getMaxChannelCount() const
{
    return 1;
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
    if (channels < getMinChannelCount() ||
        channels > getMaxChannelCount()) return false;

    m_channels = channels;
    m_blockSize = blockSize;
    m_stepSize = stepSize;
    m_frameCount = 0;
    m_spectralData.clear();
    m_numBins = blockSize / 2; // We store bins 1 to N/2 (skipping DC)

    // Initialize FFT
    if (m_fft) delete m_fft;
    m_fft = new Vamp::FFTReal(blockSize);
    m_fftOut.resize(blockSize + 2);
    
    // Reserve input buffer
    m_inputBuffer.reserve(blockSize * m_batchSize);
    
    // Reserve spectral data
    m_spectralData.reserve(10000 * m_numBins);

    // Pre-compute Hamming window
    m_window.resize(m_blockSize);
    for (size_t i = 0; i < m_blockSize; ++i) {
        m_window[i] = 0.54 - 0.46 * std::cos(2.0 * M_PI * i / (m_blockSize - 1));
    }

    return true;
}

void
ACIaccPlugin::reset()
{
    m_spectralData.clear();
    m_inputBuffer.clear();
    m_frameCount = 0;
}

Vamp::Plugin::FeatureSet
ACIaccPlugin::process(const float *const *inputBuffers, Vamp::RealTime timestamp)
{
    FeatureSet fs;
    
    // Apply window and append to input buffer
    for (size_t i = 0; i < m_blockSize; ++i) {
        m_inputBuffer.push_back(inputBuffers[0][i] * m_window[i]);
    }
    
    // Check if we have enough frames for a batch
    size_t currentFrames = m_inputBuffer.size() / m_blockSize;
    if (currentFrames >= m_batchSize) {
        processBatch(currentFrames);
        m_inputBuffer.clear();
    }
    
    return fs;
}

void ACIaccPlugin::processBatch(size_t numFrames)
{
    if (numFrames == 0) return;

    size_t blockSize = m_blockSize;
    
    // Ensure capacity
    size_t requiredSize = m_spectralData.size() + numFrames * m_numBins;
    if (m_spectralData.capacity() < requiredSize) {
        m_spectralData.reserve(std::max(requiredSize, m_spectralData.capacity() * 2));
    }

    for (size_t frame = 0; frame < numFrames; ++frame) {
        double* frameData = &m_inputBuffer[frame * blockSize];
        
        m_fft->forward(frameData, m_fftOut.data());
        
        // Store magnitude
        for (size_t i = 1; i <= blockSize / 2; ++i) {
            double real = m_fftOut[2 * i];
            double imag = m_fftOut[2 * i + 1];
            float magnitude = static_cast<float>(std::sqrt(real * real + imag * imag));
            m_spectralData.push_back(magnitude);
        }
        
        m_frameCount++;
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
    
    if (m_spectralData.empty()) return fs;

    // Calculate number of frames
    size_t numFrames = m_spectralData.size() / m_numBins;
    
    // Determine frequency range
    size_t minBinIndex = 0;
    size_t maxBinIndex = m_numBins - 1;
    
    if (m_minFreq > 0 || m_maxFreq > 0) {
        float binResolution = (m_inputSampleRate / 2.0f) / (m_blockSize / 2);
        
        if (m_minFreq > 0) {
            minBinIndex = static_cast<size_t>(m_minFreq * 1000.0f / binResolution);
        }
        if (m_maxFreq > 0) {
            maxBinIndex = static_cast<size_t>(m_maxFreq * 1000.0f / binResolution);
            if (maxBinIndex >= m_numBins) {
                maxBinIndex = m_numBins - 1;
            }
        }
    }
    
    // Calculate frames per cluster
    // m_stepSize is samples per frame (hop size)
    // m_inputSampleRate is samples per second
    // frames per second = m_inputSampleRate / m_stepSize
    float framesPerSecond = m_inputSampleRate / m_stepSize;
    size_t framesPerCluster = static_cast<size_t>(std::round(m_clusterSize * framesPerSecond));
    
    if (framesPerCluster < 1) framesPerCluster = 1;
    
    size_t numClusters = (numFrames + framesPerCluster - 1) / framesPerCluster; // Ceiling division
    
    float totalACI = 0.0f;
    
    // Iterate over clusters
    for (size_t c = 0; c < numClusters; ++c) {
        size_t startFrame = c * framesPerCluster;
        size_t endFrame = std::min(startFrame + framesPerCluster, numFrames);
        
        if (startFrame >= endFrame) continue;
        
        size_t clusterFrames = endFrame - startFrame;
        if (clusterFrames < 2) continue; // Need at least 2 frames for difference

        // Accumulators for this cluster
        size_t numActiveBins = maxBinIndex - minBinIndex + 1;
        std::vector<float> sumIntensity(numActiveBins, 0.0f);
        std::vector<float> sumAbsDiff(numActiveBins, 0.0f);
        
        // Pass 1: Calculate Sum Intensity and Sum Abs Diff for this cluster
        // Optimized loop: Time (outer) -> Frequency (inner)
        
        // First frame of the cluster
        size_t firstFrameIdx = startFrame * m_numBins;
        for (size_t b = 0; b < numActiveBins; ++b) {
            sumIntensity[b] += m_spectralData[firstFrameIdx + minBinIndex + b];
        }
        
        // Subsequent frames
        for (size_t t = startFrame; t < endFrame - 1; ++t) {
            size_t currentFrameIdx = t * m_numBins;
            size_t nextFrameIdx = (t + 1) * m_numBins;
            
            for (size_t b = 0; b < numActiveBins; ++b) {
                float valCurrent = m_spectralData[currentFrameIdx + minBinIndex + b];
                float valNext = m_spectralData[nextFrameIdx + minBinIndex + b];
                
                sumIntensity[b] += valNext;
                sumAbsDiff[b] += std::abs(valNext - valCurrent);
            }
        }
        
        // Calculate ACI for this cluster
        float clusterAci = 0.0f;
        for (size_t b = 0; b < numActiveBins; ++b) {
            if (sumIntensity[b] > 0) {
                clusterAci += sumAbsDiff[b] / sumIntensity[b];
            }
        }
        
        totalACI += clusterAci;
    }
    
    Feature f;
    f.hasTimestamp = false;
    f.values.push_back(totalACI);
    
    fs[0].push_back(f);
    
    return fs;
}
