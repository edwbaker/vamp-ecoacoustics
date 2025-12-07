// ACIPlugin.cpp - Acoustic Complexity Index
// Uses Vamp SDK FFT implementation

#include "ACIPlugin.h"
#include <cmath>
#include <algorithm>
#include <numeric>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

ACIPlugin::ACIPlugin(float inputSampleRate) :
    Plugin(inputSampleRate),
    m_minFreq(0),
    m_maxFreq(0),  // 0 means use full frequency range
    m_nbWindows(1),
    m_blockSize(0),
    m_stepSize(0),
    m_channels(0),
    m_fft(0),
    m_frameCount(0)
{
}

ACIPlugin::~ACIPlugin()
{
    delete m_fft;
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
    return "Calculates the Acoustic Complexity Index (ACI) using Vamp SDK FFT.";
}

string
ACIPlugin::getMaker() const
{
    return "ReVAMP";
}

int
ACIPlugin::getPluginVersion() const
{
    return 1;
}

string
ACIPlugin::getCopyright() const
{
    return "MIT License";
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
    return 1;
}

Vamp::Plugin::ParameterList
ACIPlugin::getParameterDescriptors() const
{
    ParameterList list;

    ParameterDescriptor d;
    
    d.identifier = "minfreq";
    d.name = "Minimum Frequency";
    d.description = "Lower frequency limit in kHz (0 = use full range)";
    d.unit = "kHz";
    d.minValue = 0;
    d.maxValue = m_inputSampleRate / 2000.0f;  // Nyquist in kHz
    d.defaultValue = 0;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "maxfreq";
    d.name = "Maximum Frequency";
    d.description = "Upper frequency limit in kHz (0 = use full range)";
    d.unit = "kHz";
    d.minValue = 0;
    d.maxValue = m_inputSampleRate / 2000.0f;
    d.defaultValue = 0;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "nbwindows";
    d.name = "Number of Windows";
    d.description = "Number of temporal windows to divide the recording";
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
    if (identifier == "minfreq") {
        return m_minFreq;
    }
    if (identifier == "maxfreq") {
        return m_maxFreq;
    }
    if (identifier == "nbwindows") {
        return static_cast<float>(m_nbWindows);
    }
    return 0;
}

void
ACIPlugin::setParameter(string identifier, float value)
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
    m_spectralData.clear();

    delete m_fft;
    m_fft = new Vamp::FFTReal(blockSize);
    m_fftOutput.resize(blockSize + 2);

    return true;
}

void
ACIPlugin::reset()
{
    m_spectralData.clear();
    m_frameCount = 0;
}

Vamp::Plugin::FeatureSet
ACIPlugin::process(const float *const *inputBuffers, Vamp::RealTime timestamp)
{
    FeatureSet fs;
    
    // Create Hamming window
    std::vector<double> window(m_blockSize);
    for (size_t i = 0; i < m_blockSize; ++i) {
        window[i] = 0.54 - 0.46 * std::cos(2.0 * M_PI * i / (m_blockSize - 1));
    }
    
    // Apply window to input
    std::vector<double> windowed(m_blockSize);
    for (size_t i = 0; i < m_blockSize; ++i) {
        windowed[i] = inputBuffers[0][i] * window[i];
    }
    
    // Compute FFT using Vamp SDK
    m_fft->forward(windowed.data(), m_fftOutput.data());
    
    // Compute magnitudes - store bins 1-256 (skip DC at bin 0)
    std::vector<float> spectrum;
    spectrum.reserve(m_blockSize / 2);

    for (size_t i = 1; i <= m_blockSize / 2; ++i) {
        double real = m_fftOutput[2 * i];
        double imag = m_fftOutput[2 * i + 1];
        double magnitude = std::sqrt(real * real + imag * imag);
        spectrum.push_back(static_cast<float>(magnitude));
    }
    
    m_spectralData.push_back(spectrum);
    m_frameCount++;
    
    return fs;
}

Vamp::Plugin::FeatureSet
ACIPlugin::getRemainingFeatures()
{
    FeatureSet fs;
    
    if (m_spectralData.empty()) return fs;

    // Handle padding frames (same logic as ACIPlugin)
    size_t paddingFrames = 0;
    if (m_stepSize > 0) {
        paddingFrames = m_blockSize / m_stepSize;
    }
    
    size_t numFrames = m_spectralData.size();
    if (numFrames > paddingFrames) {
        numFrames -= paddingFrames;
    } else {
        numFrames = 0;
    }

    // If overlapping, we might need to append a zero frame to match seewave's tail handling
    // (Logic copied from ACIPlugin fix)
    if (m_stepSize < m_blockSize && numFrames > 0) {
        // Append a frame of zeros
        std::vector<float> zeroFrame(m_spectralData[0].size(), 0.0f);
        m_spectralData.insert(m_spectralData.begin() + numFrames, zeroFrame);
        numFrames++;
    }

    // Normalize by global maximum
    float globalMax = 0.0f;
    for (size_t frame = 0; frame < numFrames; ++frame) {
        for (size_t bin = 0; bin < m_spectralData[frame].size(); ++bin) {
            if (m_spectralData[frame][bin] > globalMax) {
                globalMax = m_spectralData[frame][bin];
            }
        }
    }
    
    if (globalMax > 0) {
        for (size_t frame = 0; frame < numFrames; ++frame) {
            for (size_t bin = 0; bin < m_spectralData[frame].size(); ++bin) {
                m_spectralData[frame][bin] /= globalMax;
            }
        }
    }
    
    // Apply frequency limits
    size_t minBin = 0;
    size_t maxBin = m_blockSize / 2 - 1;
    
    if (m_minFreq > 0 || m_maxFreq > 0) {
        float binResolution = (m_inputSampleRate / 2.0f) / (m_blockSize / 2);
        
        if (m_minFreq > 0) {
            minBin = static_cast<size_t>(m_minFreq * 1000.0f / binResolution);
        }
        if (m_maxFreq > 0) {
            maxBin = static_cast<size_t>(m_maxFreq * 1000.0f / binResolution);
            if (maxBin >= m_spectralData[0].size()) {
                maxBin = m_spectralData[0].size() - 1;
            }
        }
    }
    
    // Calculate ACI
    std::vector<float> acis(m_nbWindows, 0.0f);
    
    for (int j = 0; j < m_nbWindows; ++j) {
        size_t l = numFrames;
        size_t startFrame = static_cast<size_t>(std::floor(l / m_nbWindows * j));
        size_t endFrame = static_cast<size_t>(std::floor(l / m_nbWindows * (j + 1)));
        
        if (startFrame >= endFrame) continue;
        
        for (size_t freqBin = minBin; freqBin <= maxBin; ++freqBin) {
            std::vector<float> timeSeries;
            for (size_t t = startFrame; t < endFrame; ++t) {
                timeSeries.push_back(m_spectralData[t][freqBin]);
            }
            
            if (timeSeries.size() < 2) continue;
            
            float sumIntensity = std::accumulate(timeSeries.begin(), timeSeries.end(), 0.0f);
            
            if (sumIntensity <= 0) continue;
            
            float binAci = 0.0f;
            for (size_t t = 0; t < timeSeries.size() - 1; ++t) {
                float diff = timeSeries[t + 1] - timeSeries[t];
                float normalizedDiff = diff / sumIntensity;
                binAci += std::abs(normalizedDiff);
            }
            
            acis[j] += binAci;
        }
    }
    
    float totalACI = std::accumulate(acis.begin(), acis.end(), 0.0f);
    
    Feature f;
    f.hasTimestamp = false;
    f.values.push_back(totalACI);
    
    fs[0].push_back(f);
    
    return fs;
}
