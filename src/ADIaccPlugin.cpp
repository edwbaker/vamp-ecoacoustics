// ADIaccPlugin.cpp - Acoustic Diversity Index (Accumulated)
// Matches soundecology::acoustic_diversity

#include "ADIaccPlugin.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <iomanip>
#include <map>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef ADI_DEBUG
#define ADI_DEBUG 0
#endif

ADIaccPlugin::ADIaccPlugin(float inputSampleRate) :
    Plugin(inputSampleRate),
    m_channels(0),
    m_blockSize(0),
    m_stepSize(0),
    m_numBins(0),
    m_fft(nullptr),
    m_batchSize(256),
    m_windowType(Hanning),
    m_minFreq(0.0f),
    m_maxFreq(10000.0f),
    m_dbThreshold(-50.0f),
    m_freqStep(1000.0f)
{
}

ADIaccPlugin::~ADIaccPlugin()
{
    if (m_fft) {
        delete m_fft;
        m_fft = nullptr;
    }
}

Vamp::Plugin::InputDomain
ADIaccPlugin::getInputDomain() const
{
    return TimeDomain;
}

size_t
ADIaccPlugin::getMinChannelCount() const
{
    return 1;
}

size_t
ADIaccPlugin::getMaxChannelCount() const
{
    return 2;
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
    return "Ecoacoustic-Vamp-Plugins";
}

int
ADIaccPlugin::getPluginVersion() const
{
    return 1;
}

string
ADIaccPlugin::getCopyright() const
{
    return "(C) Ed Baker 2025. Licensed under GPL (>=2).";
}

Vamp::Plugin::ParameterList
ADIaccPlugin::getParameterDescriptors() const
{
    ParameterList list;

    ParameterDescriptor d;
    d.identifier = "maxFreq";
    d.name = "Maximum Frequency";
    d.description = "Maximum frequency (kHz) to include in calculation";
    d.unit = "kHz";
    d.minValue = 0;
    d.maxValue = m_inputSampleRate / 2000.0f;
    d.defaultValue = 10.0f;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "dbThreshold";
    d.name = "dB Threshold";
    d.description = "Threshold in dBFS";
    d.unit = "dB";
    d.minValue = -120.0f;
    d.maxValue = 0.0f;
    d.defaultValue = -50.0f;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "freqStep";
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
    if (identifier == "maxFreq") return m_maxFreq;
    if (identifier == "dbThreshold") return m_dbThreshold;
    if (identifier == "freqStep") return m_freqStep;
    return 0;
}

void
ADIaccPlugin::setParameter(string identifier, float value)
{
    if (identifier == "maxFreq") {
        m_maxFreq = value;
    } else if (identifier == "dbThreshold") {
        m_dbThreshold = value;
    } else if (identifier == "freqStep") {
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
    (void)name;
}

Vamp::Plugin::OutputList
ADIaccPlugin::getOutputDescriptors() const
{
    OutputList list;

    OutputDescriptor d;
    d.identifier = "adi";
    d.name = "ADI";
    d.description = "Acoustic Diversity Index (Left channel or mono)";
    d.unit = "";
    d.hasFixedBinCount = true;
    d.binCount = 1;
    d.hasKnownExtents = false;
    d.isQuantized = false;
    d.sampleType = OutputDescriptor::VariableSampleRate;
    d.hasDuration = false;
    list.push_back(d);

    // Add right channel output for stereo files
    OutputDescriptor d2;
    d2.identifier = "adi_right";
    d2.name = "ADI (Right)";
    d2.description = "Acoustic Diversity Index (Right channel)";
    d2.unit = "";
    d2.hasFixedBinCount = true;
    d2.binCount = 1;
    d2.hasKnownExtents = false;
    d2.isQuantized = false;
    d2.sampleType = OutputDescriptor::VariableSampleRate;
    d2.hasDuration = false;
    list.push_back(d2);

    return list;
}

size_t
ADIaccPlugin::getPreferredBlockSize() const
{
    size_t wl = static_cast<size_t>(m_inputSampleRate / 10.0f);
    if (wl % 2 == 1) wl++;
    return wl;
}

size_t
ADIaccPlugin::getPreferredStepSize() const
{
    return getPreferredBlockSize();
}

bool ADIaccPlugin::initialise(size_t channels, size_t stepSize, size_t blockSize)
{
    if (channels < getMinChannelCount() ||
        channels > getMaxChannelCount()) return false;

    m_channels = channels;
    m_blockSize = blockSize;
    m_stepSize = stepSize;
    m_numBins = blockSize / 2;

    // Initialize per-channel vectors
    m_frameCount_ch.assign(m_channels, 0);
    m_globalMax_ch.assign(m_channels, 0.0f);
    m_sampleCount_ch.assign(m_channels, 0);
    m_inputBuffer_ch.resize(m_channels);
    m_spectralData_ch.resize(m_channels);
    m_bandHistograms_ch.resize(m_channels);
    
    for (size_t ch = 0; ch < m_channels; ++ch) {
        m_inputBuffer_ch[ch].reserve(blockSize * m_batchSize);
        m_spectralData_ch[ch].clear();
    }

    // Initialize FFT
    if (m_fft) delete m_fft;
    m_fft = new Vamp::FFTReal(static_cast<unsigned int>(blockSize));
    m_fftOut.resize(blockSize + 2);
    m_fftInput.resize(blockSize);
    
    m_window.resize(m_blockSize);
    double denom = static_cast<double>(m_windowType == Hamming ? m_blockSize - 1 : m_blockSize);
    if (m_windowType == Hamming) {
        for (size_t i = 0; i < m_blockSize; ++i) {
            m_window[i] = static_cast<float>(0.54 - 0.46 * std::cos(2.0 * M_PI * i / denom));
        }
    } else {
        for (size_t i = 0; i < m_blockSize; ++i) {
            m_window[i] = static_cast<float>(0.5 * (1.0 - std::cos(2.0 * M_PI * i / denom)));
        }
    }
    
    // Pre-compute bin-to-bands mapping (each bin can belong to multiple bands)
    float maxFreqHz = m_maxFreq > 0 ? m_maxFreq : m_inputSampleRate / 2.0f;
    float binResolution = (m_inputSampleRate / 2.0f) / (m_blockSize / 2);
    int numBands = static_cast<int>(std::ceil(maxFreqHz / m_freqStep));
    m_binToBands.resize(m_numBins);
    for (size_t b = 0; b < m_numBins; ++b) {
        m_binToBands[b].clear();
        int rRow = static_cast<int>(b) + 1;  // R 1-indexed
        for (int band = 0; band < numBands; ++band) {
            float bandMinHz = band * m_freqStep;
            float bandMaxHz = (band + 1) * m_freqStep;
            int miny = static_cast<int>(std::round(bandMinHz / binResolution));
            if (miny == 0) miny = 1;
            int maxy = static_cast<int>(std::round(bandMaxHz / binResolution));
            if (rRow >= miny && rRow <= maxy) {
                m_binToBands[b].push_back(band);
            }
        }
    }

    if (ADI_DEBUG) {
        std::cerr << "Window type: " << (m_windowType == Hamming ? "Hamming" : "Hanning") << std::endl;
        std::cerr << "Window[0]: " << m_window[0] << ", Window[N/2]: " << m_window[m_blockSize/2] << std::endl;
    }

    m_bandsInitialized = false;
    m_bandStartBins.clear();
    m_bandEndBins.clear();

    return true;
}

void
ADIaccPlugin::reset()
{
    for (size_t ch = 0; ch < m_channels; ++ch) {
        m_spectralData_ch[ch].clear();
        m_inputBuffer_ch[ch].clear();
        m_frameCount_ch[ch] = 0;
        m_globalMax_ch[ch] = 0.0f;
        m_sampleCount_ch[ch] = 0;
        m_bandHistograms_ch[ch].clear();
    }
    m_bandsInitialized = false;
}

Vamp::Plugin::FeatureSet
ADIaccPlugin::process(const float *const *inputBuffers, Vamp::RealTime timestamp)
{
    FeatureSet fs;
    (void)timestamp;
    
    const float* win = m_window.data();
    
    // Process each channel independently
    for (size_t ch = 0; ch < m_channels; ++ch) {
        // Track actual samples received
        m_sampleCount_ch[ch] += m_blockSize;
        
        size_t oldSize = m_inputBuffer_ch[ch].size();
        m_inputBuffer_ch[ch].resize(oldSize + m_blockSize);
        double* dst = m_inputBuffer_ch[ch].data() + oldSize;
        const float* src = inputBuffers[ch];
        for (size_t i = 0; i < m_blockSize; ++i) {
            dst[i] = static_cast<double>(src[i]) * win[i];
        }
        
        // Check if we have enough frames for a batch
        size_t currentFrames = m_inputBuffer_ch[ch].size() / m_blockSize;
        if (currentFrames >= m_batchSize) {
            processBatch(ch, currentFrames);
            m_inputBuffer_ch[ch].clear();
        }
    }
    
    return fs;
}

void ADIaccPlugin::processBatch(size_t channel, size_t numFrames)
{
    if (numFrames == 0) return;

    // Initialize bands if needed (same for all channels)
    if (!m_bandsInitialized) {
        float maxFreqHz = m_maxFreq;
        if (maxFreqHz <= 0) maxFreqHz = m_inputSampleRate / 2.0f;
        
        int numBands = static_cast<int>(std::ceil(maxFreqHz / m_freqStep));
        
        // Initialize histograms for all channels
        for (size_t ch = 0; ch < m_channels; ++ch) {
            m_bandHistograms_ch[ch].assign(numBands, std::vector<int>(30000, 0));
        }
        
        m_bandStartBins.assign(numBands, -1);
        m_bandEndBins.assign(numBands, -1);
        
        float binResolution = (m_inputSampleRate / 2.0f) / (m_blockSize / 2);
        
        for (size_t b = 0; b < m_numBins; ++b) {
            float freq = b * binResolution;
            if (freq <= maxFreqHz) {
                int bandIndex = static_cast<int>(std::floor(freq / m_freqStep));
                if (bandIndex >= 0 && bandIndex < numBands) {
                    if (m_bandStartBins[bandIndex] == -1) {
                        m_bandStartBins[bandIndex] = static_cast<int>(b);
                    }
                    m_bandEndBins[bandIndex] = static_cast<int>(b) + 1;
                }
            }
        }
        
        m_bandsInitialized = true;
    }

    size_t blockSize = m_blockSize;
    size_t numBands = m_bandHistograms_ch[channel].size();
    
    float binResolution = (m_inputSampleRate / 2.0f) / (m_blockSize / 2);
    float maxFreqHz = m_maxFreq;
    if (maxFreqHz <= 0) maxFreqHz = m_inputSampleRate / 2.0f;
    
    float& globalMax = m_globalMax_ch[channel];
    auto& histograms = m_bandHistograms_ch[channel];
    const double* fftOut = m_fftOut.data();

    for (size_t frame = 0; frame < numFrames; ++frame) {
        const double* frameData = m_inputBuffer_ch[channel].data() + frame * blockSize;
        
        // Check if this frame is all zeros (zero-padded by host)
        double frameEnergy = 0.0;
        for (size_t i = 0; i < blockSize; ++i) {
            frameEnergy += frameData[i] * frameData[i];
        }
        
        // Skip zero-padded frames
        if (frameEnergy < 1e-20) {
            continue;
        }
        
        m_fft->forward(frameData, m_fftOut.data());

        // Process bins 0 to N/2-1 (excluding Nyquist) to match seewave::spectro
        for (size_t b = 0; b < m_numBins; ++b) {
            double real = fftOut[2 * b];
            double imag = fftOut[2 * b + 1];
            double magnitude = std::sqrt(real * real + imag * imag);
            
            if (static_cast<float>(magnitude) > globalMax) globalMax = static_cast<float>(magnitude);
            
            double db = -200.0;
            if (magnitude > 1e-10) {
                db = 20.0 * std::log10(magnitude);
            }
            
            int histBin = static_cast<int>((db + 200.0) * 100.0);
            if (histBin < 0) histBin = 0;
            if (histBin >= 30000) histBin = 29999;
            
            // Use pre-computed bin-to-bands mapping instead of checking every band
            const auto& bands = m_binToBands[b];
            for (int band : bands) {
                histograms[band][histBin]++;
            }
        }
        
        m_frameCount_ch[channel]++;
    }
}

Vamp::Plugin::FeatureSet
ADIaccPlugin::getRemainingFeatures()
{
    FeatureSet fs;
    
    // Process remaining frames for each channel
    for (size_t ch = 0; ch < m_channels; ++ch) {
        if (!m_inputBuffer_ch[ch].empty()) {
            size_t remainingFrames = m_inputBuffer_ch[ch].size() / m_blockSize;
            if (remainingFrames > 0) {
                processBatch(ch, remainingFrames);
            }
            m_inputBuffer_ch[ch].clear();
        }
    }
    
    // Calculate ADI for each channel
    for (size_t ch = 0; ch < m_channels; ++ch) {
        if (m_frameCount_ch[ch] == 0) {
            Feature f;
            f.hasTimestamp = true;
            f.timestamp = Vamp::RealTime::frame2RealTime(0, static_cast<unsigned int>(m_inputSampleRate));
            f.values.push_back(0.0f);
            fs[static_cast<int>(ch)].push_back(f);
            continue;
        }

        // Calculate Threshold relative to Global Max for this channel
        double globalMaxDB = -200.0;
        if (m_globalMax_ch[ch] > 1e-10) {
            globalMaxDB = 20.0 * std::log10(m_globalMax_ch[ch]);
        }

        double thresholdDB = globalMaxDB + m_dbThreshold;
        
        int thresholdBin = static_cast<int>(std::floor((thresholdDB + 200.0) * 100.0));
        if (thresholdBin >= 30000) thresholdBin = 29999;

        
        // Calculate band proportions
        size_t numBands = m_bandHistograms_ch[ch].size();
        std::vector<double> bandProportions(numBands, 0.0);
        double sumProportions = 0.0;
        
        for (size_t i = 0; i < numBands; ++i) {
            const int* histData = m_bandHistograms_ch[ch][i].data();
            size_t histSize = m_bandHistograms_ch[ch][i].size();
            long count = 0;
            long totalCells = 0;
            
            // Sum all bins and count above threshold in single pass
            for (size_t bin = 0; bin < histSize; ++bin) {
                int val = histData[bin];
                totalCells += val;
                if (static_cast<int>(bin) > thresholdBin) {
                    count += val;
                }
            }
            
            if (totalCells > 0) {
                bandProportions[i] = static_cast<double>(count) / static_cast<double>(totalCells);
            } else {
                bandProportions[i] = 0.0;
            }
            
            sumProportions += bandProportions[i];
        }
        
        // Calculate Shannon Index
        double adi = 0.0;
        if (sumProportions > 0) {
            for (size_t i = 0; i < numBands; ++i) {
                double p = bandProportions[i] / sumProportions;
                if (p > 0) {
                    double term = -p * std::log(p);
                    adi += term;
                }
            }
        }
        
        Feature f;
        f.hasTimestamp = true;
        f.timestamp = Vamp::RealTime::frame2RealTime(0, static_cast<unsigned int>(m_inputSampleRate));
        f.values.push_back(static_cast<float>(adi));
        
        fs[static_cast<int>(ch)].push_back(f);
    }
    
    return fs;
}
