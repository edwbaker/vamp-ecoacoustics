// ACIaccPlugin.cpp - Acoustic Complexity Index (Accumulated)
// Matches soundecology::acoustic_complexity

#include "ACIaccPlugin.h"
#include <cmath>
#include <algorithm>
#include <numeric>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

ACIaccPlugin::ACIaccPlugin(float inputSampleRate) :
    Plugin(inputSampleRate),
    m_channels(0),
    m_blockSize(0),
    m_stepSize(0),
    m_numBins(0),
    m_frameCount(0),
    m_minFreq(0),
    m_maxFreq(0),
    m_fft(nullptr),
    m_windowType(Hamming),
    m_clusterSize(5.0f),
    m_framesPerCluster(0)
{
}

ACIaccPlugin::~ACIaccPlugin()
{
    if (m_fft) {
        delete m_fft;
        m_fft = nullptr;
    }
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
    return 2;
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
    return "Calculates ACI matching soundecology::acoustic_complexity().";
}

string
ACIaccPlugin::getMaker() const
{
    return "Ecoacoustic-Vamp-Plugins";
}

int
ACIaccPlugin::getPluginVersion() const
{
    return 1;
}

string
ACIaccPlugin::getCopyright() const
{
    return "(C) Ed Baker 2025. Licensed under GPL (>=2).";
}

Vamp::Plugin::ParameterList
ACIaccPlugin::getParameterDescriptors() const
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
    d.defaultValue = 0; // 0 means max
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "clusterSize";
    d.name = "Cluster Size";
    d.description = "Size of the cluster in seconds (j)";
    d.unit = "s";
    d.minValue = 0.1f;
    d.maxValue = 600.0f;
    d.defaultValue = 5.0f;
    d.isQuantized = false;
    list.push_back(d);

    return list;
}

float
ACIaccPlugin::getParameter(string identifier) const
{
    if (identifier == "minFreq") return m_minFreq;
    if (identifier == "maxFreq") return m_maxFreq;
    if (identifier == "clusterSize") return m_clusterSize;
    return 0;
}

void
ACIaccPlugin::setParameter(string identifier, float value)
{
    if (identifier == "minFreq") {
        m_minFreq = value;
    } else if (identifier == "maxFreq") {
        m_maxFreq = value;
    } else if (identifier == "clusterSize") {
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
    (void)name;
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
    m_numBins = blockSize / 2;

    // Initialize FFT
    if (m_fft) delete m_fft;
    m_fft = new Vamp::FFTReal(static_cast<unsigned int>(blockSize));
    m_fftInput.resize(blockSize);
    m_fftOut.resize(blockSize + 2);
    
    // Initialize per-channel storage
    m_spectralData_ch.resize(m_channels);
    m_totalSamplesReceived.resize(m_channels, 0);
    m_lastBlockValidSamples.resize(m_channels, 0);
    m_frameCount_ch.resize(m_channels, 0);
    m_spectralWriteIdx_ch.resize(m_channels, 0);
    
    // Calculate frames per cluster for initial reserve
    float framesPerSecond = m_inputSampleRate / m_stepSize;
    m_framesPerCluster = static_cast<size_t>(std::floor(m_clusterSize * framesPerSecond));
    if (m_framesPerCluster < 1) m_framesPerCluster = 1;
    
    for (size_t ch = 0; ch < m_channels; ++ch) {
        m_spectralData_ch[ch].reserve(4096 * m_numBins);
    }
    
    // Pre-compute window (float for memory efficiency)
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

    return true;
}

void
ACIaccPlugin::reset()
{
    for (size_t ch = 0; ch < m_channels; ++ch) {
        m_spectralData_ch[ch].clear();
        m_totalSamplesReceived[ch] = 0;
        m_lastBlockValidSamples[ch] = 0;
        m_frameCount_ch[ch] = 0;
        m_spectralWriteIdx_ch[ch] = 0;
    }
    m_frameCount = 0;
}

void ACIaccPlugin::processFrame(size_t channel, const float* input)
{
    // Apply window and copy to FFT input buffer
    const float* win = m_window.data();
    double* fftIn = m_fftInput.data();
    for (size_t i = 0; i < m_blockSize; ++i) {
        fftIn[i] = static_cast<double>(input[i]) * win[i];
    }
    
    // Perform FFT
    m_fft->forward(fftIn, m_fftOut.data());
    
    // Store magnitude squared (defer sqrt to getRemainingFeatures)
    size_t requiredSize = m_spectralWriteIdx_ch[channel] + m_numBins;
    if (m_spectralData_ch[channel].size() < requiredSize) {
        m_spectralData_ch[channel].resize(requiredSize);
    }
    
    float* spectralPtr = m_spectralData_ch[channel].data() + m_spectralWriteIdx_ch[channel];
    const double* fftOut = m_fftOut.data();
    for (size_t i = 0; i < m_numBins; ++i) {
        double real = fftOut[2 * i];
        double imag = fftOut[2 * i + 1];
        spectralPtr[i] = static_cast<float>(real * real + imag * imag);
    }
    m_spectralWriteIdx_ch[channel] = requiredSize;
    m_frameCount_ch[channel]++;
}

Vamp::Plugin::FeatureSet
ACIaccPlugin::process(const float *const *inputBuffers, Vamp::RealTime timestamp)
{
    FeatureSet fs;
    (void)timestamp;
    
    for (size_t ch = 0; ch < m_channels; ++ch) {
        // Detect trailing zeros (padding) in this block
        size_t validSamples = m_blockSize;
        for (size_t i = m_blockSize; i > 0; --i) {
            if (inputBuffers[ch][i-1] != 0.0f) {
                validSamples = i;
                break;
            }
            if (i == 1) validSamples = 0;
        }
        
        m_totalSamplesReceived[ch] += m_blockSize;
        m_lastBlockValidSamples[ch] = validSamples;
        
        processFrame(ch, inputBuffers[ch]);
    }
    
    return fs;
}

Vamp::Plugin::FeatureSet
ACIaccPlugin::getRemainingFeatures()
{
    FeatureSet fs;
    
    for (size_t ch = 0; ch < m_channels; ++ch) {
        size_t totalFrames = m_spectralWriteIdx_ch[ch] / m_numBins;
        if (totalFrames == 0) {
            Feature f;
            f.hasTimestamp = false;
            f.values.push_back(0.0f);
            fs[0].push_back(f);
            continue;
        }
        
        // Calculate duration based on true sample count (accounting for padding)
        size_t paddingSamples = m_blockSize - m_lastBlockValidSamples[ch];
        size_t trueSampleCount = m_totalSamplesReceived[ch] - paddingSamples;
        double duration = static_cast<double>(trueSampleCount) / m_inputSampleRate;
        
        // Number of clusters = floor(duration / cluster_size)
        size_t numClusters = static_cast<size_t>(std::floor(duration / m_clusterSize));
        if (numClusters == 0) numClusters = 1;
        
        // Frames per cluster (I_per_j in soundecology)
        double delta_tk = duration / totalFrames;
        size_t framesPerCluster = static_cast<size_t>(std::floor(m_clusterSize / delta_tk));
        if (framesPerCluster < 1) framesPerCluster = 1;
        
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
        
        // Get raw pointer for faster access - data is magnitude squared
        const float* data = m_spectralData_ch[ch].data();
        size_t numActiveBins = maxBinIndex - minBinIndex + 1;
        
        // Process clusters - matching soundecology's exact algorithm
        float totalACI = 0.0f;
        
        // Stack allocation for accumulators
        float stackSumI[512], stackSumD[512];
        float* sumI = (numActiveBins <= 512) ? stackSumI : new float[numActiveBins];
        float* sumD = (numActiveBins <= 512) ? stackSumD : new float[numActiveBins];
        
        // Stack buffer for previous frame magnitudes (avoid recomputing sqrt)
        float stackPrevMag[512];
        float* prevMag = (numActiveBins <= 512) ? stackPrevMag : new float[numActiveBins];
        
        for (size_t cluster = 0; cluster < numClusters; ++cluster) {
            size_t minFrame = cluster * framesPerCluster;
            size_t maxFrame = (cluster + 1) * framesPerCluster;
            
            if (maxFrame > totalFrames) break;
            
            // Reset accumulators
            std::fill(sumI, sumI + numActiveBins, 0.0f);
            std::fill(sumD, sumD + numActiveBins, 0.0f);
            
            // First frame: compute sqrt and store for reuse
            const float* firstFrameSq = data + minFrame * m_numBins + minBinIndex;
            for (size_t b = 0; b < numActiveBins; ++b) {
                float mag = std::sqrt(firstFrameSq[b]);
                sumI[b] = mag;
                prevMag[b] = mag;
            }
            
            // Process subsequent frames - compute sqrt once per bin per frame
            for (size_t t = minFrame + 1; t < maxFrame; ++t) {
                const float* frameSq = data + t * m_numBins + minBinIndex;
                
                for (size_t b = 0; b < numActiveBins; ++b) {
                    float mag = std::sqrt(frameSq[b]);
                    sumI[b] += mag;
                    sumD[b] += std::abs(mag - prevMag[b]);
                    prevMag[b] = mag;
                }
            }
            
            // ACI for this cluster
            for (size_t b = 0; b < numActiveBins; ++b) {
                if (sumI[b] > 0.0f) {
                    totalACI += sumD[b] / sumI[b];
                }
            }
        }
        
        if (numActiveBins > 512) {
            delete[] sumI;
            delete[] sumD;
            delete[] prevMag;
        }

        Feature f;
        f.hasTimestamp = false;
        f.values.push_back(totalACI);
        fs[0].push_back(f);
    }
    
    return fs;
}
