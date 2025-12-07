// ACImtPlugin.cpp - Acoustic Complexity Index (Fast FFT)

#include "ACImtPlugin.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

ACImtPlugin::ACImtPlugin(float inputSampleRate) :
    Plugin(inputSampleRate),
    m_minFreq(0),
    m_maxFreq(0),
    m_nbWindows(1),
    m_channels(0),
    m_blockSize(0),
    m_stepSize(0),
    m_frameCount(0),
    m_fft(nullptr),
    m_batchSize(256) // Process 256 frames at a time
{
}

ACImtPlugin::~ACImtPlugin()
{
    if (m_fft) {
        delete m_fft;
        m_fft = nullptr;
    }
}

string
ACImtPlugin::getIdentifier() const
{
    return "acimt";
}

string
ACImtPlugin::getName() const
{
    return "Acoustic Complexity Index (BLAS MT)";
}

string
ACImtPlugin::getDescription() const
{
    return "Calculates ACI using BLAS matrix multiplication for spectrogram generation.";
}

string
ACImtPlugin::getMaker() const
{
    return "ReVAMP";
}

int
ACImtPlugin::getPluginVersion() const
{
    return 1;
}

string
ACImtPlugin::getCopyright() const
{
    return "MIT License";
}

Vamp::Plugin::InputDomain
ACImtPlugin::getInputDomain() const
{
    return TimeDomain;
}

size_t
ACImtPlugin::getPreferredBlockSize() const
{
    return 512;
}

size_t
ACImtPlugin::getPreferredStepSize() const
{
    return 512;
}

size_t
ACImtPlugin::getMinChannelCount() const
{
    return 1;
}

size_t
ACImtPlugin::getMaxChannelCount() const
{
    return 1;
}

Vamp::Plugin::ParameterList
ACImtPlugin::getParameterDescriptors() const
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

    d.identifier = "nbwindows";
    d.name = "Number of Windows";
    d.description = "Number of time windows to divide the file into";
    d.unit = "";
    d.minValue = 1;
    d.maxValue = 100;
    d.defaultValue = 1;
    d.isQuantized = true;
    d.quantizeStep = 1;
    list.push_back(d);

    return list;
}

float
ACImtPlugin::getParameter(string identifier) const
{
    if (identifier == "minfreq") return m_minFreq;
    if (identifier == "maxfreq") return m_maxFreq;
    if (identifier == "nbwindows") return m_nbWindows;
    return 0;
}

void
ACImtPlugin::setParameter(string identifier, float value)
{
    if (identifier == "minfreq") {
        m_minFreq = value;
    } else if (identifier == "maxfreq") {
        m_maxFreq = value;
    } else if (identifier == "nbwindows") {
        m_nbWindows = static_cast<int>(value);
    }
}

Vamp::Plugin::ProgramList
ACImtPlugin::getPrograms() const
{
    ProgramList list;
    return list;
}

string
ACImtPlugin::getCurrentProgram() const
{
    return "";
}

void
ACImtPlugin::selectProgram(string name)
{
}

Vamp::Plugin::OutputList
ACImtPlugin::getOutputDescriptors() const
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
ACImtPlugin::initialise(size_t channels, size_t stepSize, size_t blockSize)
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
    
    // Reserve spectral data (heuristic: assume 1 minute of audio at 44.1kHz with 512 hop)
    // ~5000 frames * 256 bins = 1.2M floats = 5MB. Safe to reserve some.
    m_spectralData.reserve(10000 * m_numBins);

    // Pre-compute Hamming window
    m_window.resize(m_blockSize);
    for (size_t i = 0; i < m_blockSize; ++i) {
        m_window[i] = 0.54 - 0.46 * std::cos(2.0 * M_PI * i / (m_blockSize - 1));
    }

    return true;
}

void
ACImtPlugin::reset()
{
    m_spectralData.clear();
    m_inputBuffer.clear();
    m_frameCount = 0;
}

Vamp::Plugin::FeatureSet
ACImtPlugin::process(const float *const *inputBuffers, Vamp::RealTime timestamp)
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

void ACImtPlugin::processBatch(size_t numFrames)
{
    if (numFrames == 0) return;

    size_t blockSize = m_blockSize;
    // We store bins 1 to blockSize/2. Total bins stored = blockSize/2.
    
    // Ensure capacity
    size_t requiredSize = m_spectralData.size() + numFrames * m_numBins;
    if (m_spectralData.capacity() < requiredSize) {
        m_spectralData.reserve(std::max(requiredSize, m_spectralData.capacity() * 2));
    }

    for (size_t frame = 0; frame < numFrames; ++frame) {
        // Get pointer to current frame in input buffer
        double* frameData = &m_inputBuffer[frame * blockSize];
        
        // Perform FFT
        m_fft->forward(frameData, m_fftOut.data());
        
        // Compute magnitudes and append to flattened vector
        // m_fftOut contains interleaved complex data: real, imag, real, imag...
        // Skip DC (i=0), start from i=1 up to Nyquist (i=blockSize/2)
        for (size_t i = 1; i <= blockSize / 2; ++i) {
            double real = m_fftOut[2 * i];
            double imag = m_fftOut[2 * i + 1];
            // Use float for storage to save memory/bandwidth
            float magnitude = static_cast<float>(std::sqrt(real * real + imag * imag));
            m_spectralData.push_back(magnitude);
        }
        
        m_frameCount++;
    }
}

Vamp::Plugin::FeatureSet
ACImtPlugin::getRemainingFeatures()
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

    // Handle padding frames
    size_t paddingFrames = 0;
    if (m_stepSize > 0) {
        paddingFrames = m_blockSize / m_stepSize;
    }
    
    // Calculate number of frames currently stored
    size_t numFrames = m_spectralData.size() / m_numBins;
    
    if (numFrames > paddingFrames) {
        numFrames -= paddingFrames;
    } else {
        numFrames = 0;
    }

    if (m_stepSize < m_blockSize && numFrames > 0) {
        // Append a frame of zeros (flattened)
        // m_numBins zeros
        m_spectralData.insert(m_spectralData.begin() + numFrames * m_numBins, m_numBins, 0.0f);
        numFrames++;
    }

    // Normalize by global maximum
    // Optimized: Single pass over the flattened vector
    // Only up to numFrames * m_numBins
    size_t totalElements = numFrames * m_numBins;
    
    float globalMax = 0.0f;
    for (size_t i = 0; i < totalElements; ++i) {
        if (m_spectralData[i] > globalMax) {
            globalMax = m_spectralData[i];
        }
    }
    
    if (globalMax > 0) {
        float invMax = 1.0f / globalMax;
        for (size_t i = 0; i < totalElements; ++i) {
            m_spectralData[i] *= invMax;
        }
    }
    
    // Apply frequency limits
    // m_spectralData is flattened: [Frame0_Bin0, Frame0_Bin1... Frame1_Bin0...]
    // m_numBins is the number of bins per frame (blockSize/2)
    
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
    
    // Calculate ACI
    // Optimized: Iterate Time (outer) -> Frequency (inner) for cache locality
    // But ACI requires summing differences per frequency bin first?
    // ACI_tot = Sum_k ( Sum_t |I_k(t+1) - I_k(t)| / Sum_t I_k(t) )
    // We need Sum_t I_k(t) (Total Intensity per bin) and Sum_t |Diff| (Total Difference per bin)
    
    std::vector<float> acis(m_nbWindows, 0.0f);
    
    for (int j = 0; j < m_nbWindows; ++j) {
        size_t l = numFrames;
        size_t startFrame = static_cast<size_t>(std::floor(l / m_nbWindows * j));
        size_t endFrame = static_cast<size_t>(std::floor(l / m_nbWindows * (j + 1)));
        
        if (startFrame >= endFrame) continue;
        
        size_t numSteps = endFrame - startFrame;
        if (numSteps < 2) continue;

        // Accumulators for each bin in the frequency range
        size_t numActiveBins = maxBinIndex - minBinIndex + 1;
        std::vector<float> sumIntensity(numActiveBins, 0.0f);
        std::vector<float> sumAbsDiff(numActiveBins, 0.0f);
        
        // Pass 1: Calculate Sum Intensity and Sum Abs Diff
        // We iterate frames (outer) to be cache friendly with the flattened vector
        
        // First frame of the window
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
                
                sumIntensity[b] += valNext; // Add next frame to intensity sum
                sumAbsDiff[b] += std::abs(valNext - valCurrent);
            }
        }
        
        // Finalize ACI for this window
        float windowAci = 0.0f;
        for (size_t b = 0; b < numActiveBins; ++b) {
            if (sumIntensity[b] > 0) {
                windowAci += sumAbsDiff[b] / sumIntensity[b];
            }
        }
        
        acis[j] = windowAci;
    }
    
    float totalACI = std::accumulate(acis.begin(), acis.end(), 0.0f);
    
    Feature f;
    f.hasTimestamp = false;
    f.values.push_back(totalACI);
    
    fs[0].push_back(f);
    
    return fs;
}
