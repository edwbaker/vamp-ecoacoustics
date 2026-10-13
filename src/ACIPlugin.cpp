/*
    Ecoacoustic Vamp Plugins
    
    High-performance implementations of acoustic indices for 
    bioacoustics and soundscape ecology analysis.
    
    (C) Ed Baker 2025. Licensed under GPL (>=2).
*/

// ACIPlugin.cpp - Acoustic Complexity Index (Fast FFT)
// Matches seewave::ACI

#include "ACIPlugin.h"
#include "RealFFT.h"
#include <cmath>
#include <algorithm>
#include <numeric>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

ACIPlugin::ACIPlugin(float inputSampleRate) :
    Plugin(inputSampleRate),
    m_nbWindows(1),
    m_channels(0),
    m_blockSize(0),
    m_stepSize(0),
    m_numBins(0),
    m_frameCount(0),
    m_globalMaxSq(0.0f),
    m_minFreq(0),
    m_maxFreq(0),
    m_fft(nullptr),
    m_hasFirstFrame(false),
    m_windowType(Hamming)
{
}

ACIPlugin::~ACIPlugin()
{
    if (m_fft) {
        delete m_fft;
        m_fft = nullptr;
    }
}

Vamp::Plugin::InputDomain
ACIPlugin::getInputDomain() const
{
    return TimeDomain;
}

size_t
ACIPlugin::getPreferredBlockSize() const
{
    return 512;
}

size_t
ACIPlugin::getPreferredStepSize() const
{
    return 512;
}

size_t
ACIPlugin::getMinChannelCount() const
{
    return 1;
}

size_t
ACIPlugin::getMaxChannelCount() const
{
    return 2;  // Accept stereo - we use channel 1 only, matching seewave
}

string
ACIPlugin::getIdentifier() const
{
    return "aci";
}

string
ACIPlugin::getName() const
{
    return "Acoustic Complexity Index";
}

string
ACIPlugin::getDescription() const
{
    return "Calculates the Acoustic Complexity Index (ACI) of a signal.";
}

string
ACIPlugin::getMaker() const
{
    return "Ecoacoustic-Vamp-Plugins";
}

int
ACIPlugin::getPluginVersion() const
{
    return 1;
}

string
ACIPlugin::getCopyright() const
{
    return "(C) Ed Baker 2025. Licensed under GPL (>=2).";
}

Vamp::Plugin::ParameterList
ACIPlugin::getParameterDescriptors() const
{
    ParameterList list;

    ParameterDescriptor d;
    d.identifier = "minFreq";
    d.name = "Minimum Frequency";
    d.description = "Minimum frequency (kHz) to include in calculation";
    d.unit = "kHz";
    d.minValue = 0;
    d.maxValue = m_inputSampleRate / 2000.0f;
    d.defaultValue = 0;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "maxFreq";
    d.name = "Maximum Frequency";
    d.description = "Maximum frequency (kHz) to include in calculation";
    d.unit = "kHz";
    d.minValue = 0;
    d.maxValue = m_inputSampleRate / 2000.0f;
    d.defaultValue = 0;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "nbWindows";
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
ACIPlugin::getParameter(string identifier) const
{
    if (identifier == "minFreq") return m_minFreq;
    if (identifier == "maxFreq") return m_maxFreq;
    if (identifier == "nbWindows") return static_cast<float>(m_nbWindows);
    return 0;
}

void
ACIPlugin::setParameter(string identifier, float value)
{
    if (identifier == "minFreq") {
        m_minFreq = value;
    } else if (identifier == "maxFreq") {
        m_maxFreq = value;
    } else if (identifier == "nbWindows") {
        m_nbWindows = static_cast<int>(value);
    }
}

Vamp::Plugin::ProgramList
ACIPlugin::getPrograms() const
{
    ProgramList list;
    return list;
}

string
ACIPlugin::getCurrentProgram() const
{
    return "";
}

void
ACIPlugin::selectProgram(string name)
{
    (void)name;
}

Vamp::Plugin::OutputList
ACIPlugin::getOutputDescriptors() const
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
ACIPlugin::initialise(size_t channels, size_t stepSize, size_t blockSize)
{
    if (channels < getMinChannelCount() ||
        channels > getMaxChannelCount()) return false;

    m_channels = channels;
    m_blockSize = blockSize;
    m_stepSize = stepSize;
    m_frameCount = 0;
    m_globalMaxSq = 0.0f;
    m_numBins = blockSize / 2 + 1;
    m_hasFirstFrame = false;
    m_spectral.init(m_numBins);

    // Initialize FFT
    if (m_fft) delete m_fft;
    m_fft = new InPlaceRealFFT(blockSize);
    
    // Pre-compute window
    m_window.resize(m_blockSize);
    double denom = static_cast<double>(m_blockSize - 1);
    if (m_windowType == Hamming) {
        for (size_t i = 0; i < m_blockSize; ++i) {
            m_window[i] = static_cast<float>(0.54 - 0.46 * std::cos(2.0 * M_PI * i / denom));
        }
    } else {
        for (size_t i = 0; i < m_blockSize; ++i) {
            m_window[i] = static_cast<float>(0.5 * (1.0 - std::cos(2.0 * M_PI * i / denom)));
        }
    }

    // Streaming mode buffers (for single window - most common case)
    m_prevFrame.resize(m_numBins, 0.0f);
    m_sumIntensity.resize(m_numBins, 0.0f);
    m_sumAbsDiff.resize(m_numBins, 0.0f);

    return true;
}

void
ACIPlugin::reset()
{
    m_spectral.clear();
    m_frameCount = 0;
    m_globalMaxSq = 0.0f;
    m_hasFirstFrame = false;
    std::fill(m_prevFrame.begin(), m_prevFrame.end(), 0.0f);
    std::fill(m_sumIntensity.begin(), m_sumIntensity.end(), 0.0f);
    std::fill(m_sumAbsDiff.begin(), m_sumAbsDiff.end(), 0.0f);
}

void ACIPlugin::processFrame(const float* input)
{
    // Apply window, writing the samples straight into the FFT's own buffer
    double* fftIn = m_fft->buffer();
    for (size_t i = 0; i < m_blockSize; ++i) {
        fftIn[i] = static_cast<double>(input[i] * m_window[i]);
    }

    // Perform FFT
    m_fft->forward();

    size_t halfBlock = m_blockSize / 2;

    // Store spectral data (magnitude squared)
    float* spectralPtr = m_spectral.appendRow();
    float globalMaxSq = m_globalMaxSq;
    m_fft->forBins(0, halfBlock + 1, [spectralPtr, &globalMaxSq](size_t i, double real, double imag) {
        float magSq = static_cast<float>(real * real + imag * imag);
        spectralPtr[i] = magSq;
        if (magSq > globalMaxSq) globalMaxSq = magSq;
    });
    m_globalMaxSq = globalMaxSq;
    m_frameCount++;
}

Vamp::Plugin::FeatureSet
ACIPlugin::process(const float *const *inputBuffers, Vamp::RealTime timestamp)
{
    FeatureSet fs;
    (void)timestamp;
    
    processFrame(inputBuffers[0]);
    
    return fs;
}

Vamp::Plugin::FeatureSet
ACIPlugin::getRemainingFeatures()
{
    FeatureSet fs;
    
    if (m_spectral.rows() == 0) return fs;

    // Calculate number of frames currently stored
    size_t numFrames = m_spectral.rows();
    
    // Handle padding frames
    size_t paddingFrames = 0;
    if (m_stepSize > 0) {
        paddingFrames = m_blockSize / m_stepSize;
    }
    
    if (numFrames > paddingFrames) {
        numFrames -= paddingFrames;
    } else {
        numFrames = 0;
    }
    
    if (numFrames < 2) return fs;

    // Apply frequency limits 
    size_t minBinIndex = 1; // Default to 1 to exclude DC
    size_t maxBinIndex = m_numBins - 1;
    
    if (m_minFreq > 0 || m_maxFreq > 0) {
        float binResolution = (m_inputSampleRate / 2.0f) / (m_blockSize / 2);
        
        if (m_minFreq > 0) {
            size_t calculatedMin = static_cast<size_t>(m_minFreq * 1000.0f / binResolution);
            if (calculatedMin > minBinIndex) {
                minBinIndex = calculatedMin;
            }
        }
        if (m_maxFreq > 0) {
            maxBinIndex = static_cast<size_t>(m_maxFreq * 1000.0f / binResolution);
            if (maxBinIndex >= m_numBins) {
                maxBinIndex = m_numBins - 1;
            }
        }
    }
    
    // The stored spectra are magnitude squared. A frame is turned into normalized magnitude, sqrt(x) * invMax with
    // invMax from the global maximum seen while processing, as it is read below, rather than in a separate pass
    // over everything first. The all-zero frame added next (frame numFrames) is used as it is.
    const bool normalize = m_globalMaxSq > 0;
    const float invMax = normalize ? 1.0f / std::sqrt(m_globalMaxSq) : 0.0f;
    const size_t normalizedFrames = numFrames;

    if (m_stepSize < m_blockSize) {
        // One extra, all-zero frame follows the last one used (it takes the place of the first padding frame)
        float* zeroFrame = (numFrames < m_spectral.rows()) ? m_spectral.row(numFrames) : m_spectral.appendRow();
        std::fill(zeroFrame, zeroFrame + m_numBins, 0.0f);
        numFrames++;
    }
    
    float totalACI = 0.0f;
    size_t numActiveBins = maxBinIndex - minBinIndex + 1;

    auto loadFrame = [&](size_t t, float* out) {
        const float* in = m_spectral.row(t) + minBinIndex;
        if (normalize && t < normalizedFrames) {
            for (size_t b = 0; b < numActiveBins; ++b) {
                out[b] = std::sqrt(in[b]) * invMax;
            }
        } else {
            for (size_t b = 0; b < numActiveBins; ++b) {
                out[b] = in[b];
            }
        }
    };
    
    for (int j = 0; j < m_nbWindows; ++j) {
        size_t l = numFrames;
        size_t startFrame = static_cast<size_t>(l / m_nbWindows * j);
        size_t endFrame = static_cast<size_t>(l / m_nbWindows * (j + 1));
        
        if (startFrame >= endFrame || (endFrame - startFrame) < 2) continue;

        // Stack allocation for accumulators
        float stackSumInt[512], stackSumDiff[512], stackFrameA[512], stackFrameB[512];
        float* sumIntensity = (numActiveBins <= 512) ? stackSumInt : new float[numActiveBins];
        float* sumAbsDiff = (numActiveBins <= 512) ? stackSumDiff : new float[numActiveBins];
        float* currentFrame = (numActiveBins <= 512) ? stackFrameA : new float[numActiveBins];
        float* nextFrame = (numActiveBins <= 512) ? stackFrameB : new float[numActiveBins];
        
        std::fill(sumIntensity, sumIntensity + numActiveBins, 0.0f);
        std::fill(sumAbsDiff, sumAbsDiff + numActiveBins, 0.0f);
        
        loadFrame(startFrame, currentFrame);
        for (size_t b = 0; b < numActiveBins; ++b) {
            sumIntensity[b] = currentFrame[b];
        }

        for (size_t t = startFrame; t < endFrame - 1; ++t) {
            loadFrame(t + 1, nextFrame);
            
            for (size_t b = 0; b < numActiveBins; ++b) {
                sumIntensity[b] += nextFrame[b];
                sumAbsDiff[b] += std::abs(nextFrame[b] - currentFrame[b]);
            }
            std::swap(currentFrame, nextFrame);
        }
        
        float windowAci = 0.0f;
        for (size_t b = 0; b < numActiveBins; ++b) {
            if (sumIntensity[b] > 0) {
                windowAci += sumAbsDiff[b] / sumIntensity[b];
            }
        }
        
        totalACI += windowAci;
        
        if (numActiveBins > 512) {
            delete[] sumIntensity;
            delete[] sumAbsDiff;
            delete[] currentFrame;
            delete[] nextFrame;
        }
    }
    
    // The spectra are not needed any more: release them now rather than with the plugin (other plugins run
    // together with this one may still have a large FFT to do)
    m_spectral.clear();
    
    Feature f;
    f.hasTimestamp = false;
    f.values.push_back(totalACI);
    
    fs[0].push_back(f);
    
    return fs;
}

