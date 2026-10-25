/*
    Ecoacoustic Vamp Plugins
    
    High-performance implementations of acoustic indices for 
    bioacoustics and soundscape ecology analysis.
    
    (C) Ed Baker 2025. Licensed under GPL (>=2).
*/

// ADIaccPlugin.cpp - Acoustic Diversity Index (Accumulated)
// Matches soundecology::acoustic_diversity

#include "ADIaccPlugin.h"
#include "RealFFT.h"
#include "DbAbove.h"
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
    m_windowType(Hanning),
    m_minFreq(0.0f),
    m_maxFreq(10000.0f),
    m_dbThreshold(-50.0f),
    m_freqStep(1000.0f),
    m_storedBins(0),
    m_numBands(0)
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
    m_spectralData_ch.resize(m_channels);
    m_magnitude_ch.resize(m_channels);
    m_maxBeforeLast_ch.assign(m_channels, 0.0f);
    m_lastPadded_ch.assign(m_channels, 0);

    for (size_t ch = 0; ch < m_channels; ++ch) {
        m_spectralData_ch[ch].clear();
    }

    // Initialize FFT
    if (m_fft) delete m_fft;
    m_fft = new InPlaceRealFFT(blockSize);

    // Symmetric window, as seewave's hanning.w and hamming.w (denominator n - 1)
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
    m_numBands = numBands;

    // Only the bins that belong to a band are kept for the counts
    m_storedBins = 0;
    for (size_t b = 0; b < m_numBins; ++b) {
        if (!m_binToBands[b].empty()) m_storedBins = b + 1;
    }
    for (size_t ch = 0; ch < m_channels; ++ch) {
        m_magnitude_ch[ch].init(m_storedBins);
    }

    if (ADI_DEBUG) {
        std::cerr << "Window type: " << (m_windowType == Hamming ? "Hamming" : "Hanning") << std::endl;
        std::cerr << "Window[0]: " << m_window[0] << ", Window[N/2]: " << m_window[m_blockSize/2] << std::endl;
    }

    return true;
}

void
ADIaccPlugin::reset()
{
    for (size_t ch = 0; ch < m_channels; ++ch) {
        m_spectralData_ch[ch].clear();
        m_frameCount_ch[ch] = 0;
        m_globalMax_ch[ch] = 0.0f;
        m_sampleCount_ch[ch] = 0;
        m_magnitude_ch[ch].clear();
        m_maxBeforeLast_ch[ch] = 0.0f;
        m_lastPadded_ch[ch] = 0;
    }
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
        
        // Window the block straight into the FFT's own buffer and transform it at once
        double* dst = m_fft->buffer();
        const float* src = inputBuffers[ch];
        for (size_t i = 0; i < m_blockSize; ++i) {
            dst[i] = static_cast<double>(src[i]) * win[i];
        }

        // A block ending in zeros is the host's padding of an incomplete last window (the rule ACIaccPlugin uses)
        processFrame(ch, src[m_blockSize - 1] == 0.0f);
    }
    
    return fs;
}

void ADIaccPlugin::processFrame(size_t channel, bool padded)
{
    size_t blockSize = m_blockSize;

    double& globalMax = m_globalMax_ch[channel];
    const double* frameData = m_fft->buffer();

    // Check if this frame is all zeros (zero-padded by host). The energy is a sum of non-negative terms, so it is
    // at least as large as any one of them: once a single term reaches the threshold the frame is not a padded one,
    // whatever the rest add up to, and the sum need not be finished. Otherwise the sum is the same as always.
    double frameEnergy = 0.0;
    for (size_t i = 0; i < blockSize; ++i) {
        double term = frameData[i] * frameData[i];
        if (term >= 1e-20) {
            frameEnergy = term;
            break;
        }
        frameEnergy += term;
    }
    
    // Skip zero-padded frames
    if (frameEnergy < 1e-20) {
        return;
    }
    
    m_fft->forward();

    m_maxBeforeLast_ch[channel] = globalMax;
    m_lastPadded_ch[channel] = padded;

    // Process bins 0 to N/2-1 (excluding Nyquist) to match seewave::spectro. The global maximum is over all of them;
    // only the bins that belong to a band are kept for the counts.
    double* magnitudes = m_magnitude_ch[channel].appendRow();
    const size_t storedBins = m_storedBins;
    m_fft->forBins(0, m_numBins, [magnitudes, storedBins, &globalMax](size_t b, double real, double imag) {
        double magnitude = std::sqrt(real * real + imag * imag);

        if (magnitude > globalMax) globalMax = magnitude;

        if (b < storedBins) magnitudes[b] = magnitude;
    });

    m_frameCount_ch[channel]++;
}

Vamp::Plugin::FeatureSet
ADIaccPlugin::getRemainingFeatures()
{
    FeatureSet fs;
    
    // Calculate ADI for each channel
    for (size_t ch = 0; ch < m_channels; ++ch) {
        // seewave::spectro uses only complete windows: leave out a last frame whose block the host padded with
        // zeros, unless it is the only frame
        size_t numFrames = m_magnitude_ch[ch].rows();
        if (m_lastPadded_ch[ch] && numFrames > 1) {
            --numFrames;
            m_globalMax_ch[ch] = m_maxBeforeLast_ch[ch];
        }

        if (numFrames == 0) {
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

        // Count the cells above the threshold in each band, as AEIaccPlugin does: a cell is above the threshold when
        // 20 * log10(magnitude) > thresholdDB (-200 dB for magnitudes up to 1e-10); DbAbove gives that answer and
        // calls log10 only for magnitudes close to the threshold (see DbAbove.h).
        size_t numBands = static_cast<size_t>(m_numBands);
        std::vector<long> aboveCount(numBands, 0);
        std::vector<long> totalCount(numBands, 0);
        const DbAbove isAbove(thresholdDB);

        for (size_t frame = 0; frame < numFrames; ++frame) {
            const double* magnitudes = m_magnitude_ch[ch].row(frame);
            for (size_t b = 0; b < m_storedBins; ++b) {
                const auto& bands = m_binToBands[b];
                if (bands.empty()) continue;
                const bool above = isAbove(magnitudes[b]);
                for (int band : bands) {
                    ++totalCount[band];
                    if (above) ++aboveCount[band];
                }
            }
        }

        m_magnitude_ch[ch].clear();

        // Calculate band proportions
        std::vector<double> bandProportions(numBands, 0.0);
        double sumProportions = 0.0;

        for (size_t i = 0; i < numBands; ++i) {
            if (totalCount[i] > 0) {
                bandProportions[i] = static_cast<double>(aboveCount[i]) / static_cast<double>(totalCount[i]);
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
