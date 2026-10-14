/*
    Ecoacoustic Vamp Plugins
    
    High-performance implementations of acoustic indices for 
    bioacoustics and soundscape ecology analysis.
    
    (C) Ed Baker 2025. Licensed under GPL (>=2).
*/

// AEIaccPlugin.cpp - Acoustic Evenness Index (Accumulated)
// Matches soundecology::acoustic_evenness

#include "AEIaccPlugin.h"
#include "RealFFT.h"
#include "DbAbove.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace {
// Declared defaults, shared by the constructor and getParameterDescriptors()
const float kDefaultMaxFreqKHz = 10.0f;   // kHz (held internally in Hz)
const float kDefaultDbThreshold = -50.0f; // dBFS
const float kDefaultFreqStep = 1000.0f;   // Hz
}

AEIaccPlugin::AEIaccPlugin(float inputSampleRate) :
    Plugin(inputSampleRate),
    m_channels(0),
    m_blockSize(0),
    m_stepSize(0),
    m_numBins(0),
    m_fft(nullptr),
    m_windowType(Hanning),
    m_minFreq(0.0f),
    m_maxFreq(kDefaultMaxFreqKHz * 1000.0f),
    m_dbThreshold(kDefaultDbThreshold),
    m_freqStep(kDefaultFreqStep)
{
}

AEIaccPlugin::~AEIaccPlugin()
{
    if (m_fft) {
        delete m_fft;
        m_fft = nullptr;
    }
}

Vamp::Plugin::InputDomain
AEIaccPlugin::getInputDomain() const
{
    return TimeDomain;
}

size_t
AEIaccPlugin::getMinChannelCount() const
{
    return 1;
}

size_t
AEIaccPlugin::getMaxChannelCount() const
{
    return 2;
}

string
AEIaccPlugin::getIdentifier() const
{
    return "aei-acc";
}

string
AEIaccPlugin::getName() const
{
    return "Acoustic Evenness Index (soundecology)";
}

string
AEIaccPlugin::getDescription() const
{
    return "Calculates the Acoustic Evenness Index (AEI) of a signal, matching soundecology implementation.";
}

string
AEIaccPlugin::getMaker() const
{
    return "Ecoacoustic-Vamp-Plugins";
}

int
AEIaccPlugin::getPluginVersion() const
{
    return 1;
}

string
AEIaccPlugin::getCopyright() const
{
    return "(C) Ed Baker 2025. Licensed under GPL (>=2).";
}

Vamp::Plugin::ParameterList
AEIaccPlugin::getParameterDescriptors() const
{
    ParameterList list;

    ParameterDescriptor d;
    d.identifier = "maxFreq";
    d.name = "Maximum Frequency";
    d.description = "Maximum frequency (kHz) to include in calculation";
    d.unit = "kHz";
    d.minValue = 0;
    d.maxValue = m_inputSampleRate / 2000.0f;
    d.defaultValue = kDefaultMaxFreqKHz;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "dbThreshold";
    d.name = "dB Threshold";
    d.description = "Threshold in dBFS";
    d.unit = "dB";
    d.minValue = -120.0f;
    d.maxValue = 0.0f;
    d.defaultValue = kDefaultDbThreshold;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "freqStep";
    d.name = "Frequency Step";
    d.description = "Size of frequency bands (Hz)";
    d.unit = "Hz";
    d.minValue = 100.0f;
    d.maxValue = 10000.0f;
    d.defaultValue = kDefaultFreqStep;
    d.isQuantized = false;
    list.push_back(d);

    return list;
}

float
AEIaccPlugin::getParameter(string identifier) const
{
    if (identifier == "maxFreq") return m_maxFreq / 1000.0f;
    if (identifier == "dbThreshold") return m_dbThreshold;
    if (identifier == "freqStep") return m_freqStep;
    return 0;
}

void
AEIaccPlugin::setParameter(string identifier, float value)
{
    if (identifier == "maxFreq") {
        m_maxFreq = value * 1000.0f;
    } else if (identifier == "dbThreshold") {
        m_dbThreshold = value;
    } else if (identifier == "freqStep") {
        m_freqStep = value;
    }
}

Vamp::Plugin::ProgramList
AEIaccPlugin::getPrograms() const
{
    ProgramList list;
    return list;
}

string
AEIaccPlugin::getCurrentProgram() const
{
    return "";
}

void
AEIaccPlugin::selectProgram(string name)
{
    (void)name;
}

Vamp::Plugin::OutputList
AEIaccPlugin::getOutputDescriptors() const
{
    OutputList list;

    OutputDescriptor d;
    d.identifier = "aei";
    d.name = "AEI";
    d.description = "Acoustic Evenness Index (Left channel or mono)";
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
    d2.identifier = "aei_right";
    d2.name = "AEI (Right)";
    d2.description = "Acoustic Evenness Index (Right channel)";
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
AEIaccPlugin::getPreferredBlockSize() const
{
    size_t size = static_cast<size_t>(m_inputSampleRate / 10.0f);
    if (size % 2 != 0) size++;
    return size;
}

size_t
AEIaccPlugin::getPreferredStepSize() const
{
    return getPreferredBlockSize();
}

bool
AEIaccPlugin::initialise(size_t channels, size_t stepSize, size_t blockSize)
{
    if (channels < getMinChannelCount() ||
        channels > getMaxChannelCount()) return false;

    m_channels = channels;
    m_blockSize = blockSize;
    m_stepSize = stepSize;
    m_numBins = blockSize / 2 + 1; // We store bins 0 to N/2 (including DC)

    // Initialize per-channel vectors
    m_frameCount_ch.assign(m_channels, 0);
    m_globalMax_ch.assign(m_channels, 0.0f);
    m_magnitude_ch.resize(m_channels);
    m_spectralData_ch.resize(m_channels);
    m_bandAboveThresholdCount_ch.resize(m_channels);
    m_bandTotalCount_ch.resize(m_channels);
    
    for (size_t ch = 0; ch < m_channels; ++ch) {
        m_magnitude_ch[ch].init(m_numBins);
        m_spectralData_ch[ch].clear();
    }

    // Initialize FFT
    if (m_fft) delete m_fft;
    m_fft = new InPlaceRealFFT(blockSize);
    
    // Pre-compute window
    m_window.resize(m_blockSize);

    for (size_t i = 0; i < m_blockSize; ++i) {
        m_window[i] = 0.5 * (1.0 - std::cos(2.0 * M_PI * i / (m_blockSize - 1))); //Apply Hanning window
    }

    m_bandsInitialized = false;
    m_bandStartBins.clear();
    m_bandEndBins.clear();

    return true;
}

void
AEIaccPlugin::reset()
{
    for (size_t ch = 0; ch < m_channels; ++ch) {
        m_spectralData_ch[ch].clear();
        m_magnitude_ch[ch].clear();
        m_frameCount_ch[ch] = 0;
        m_globalMax_ch[ch] = 0.0f;
        m_bandAboveThresholdCount_ch[ch].clear();
        m_bandTotalCount_ch[ch].clear();
    }
    m_bandsInitialized = false;
}

Vamp::Plugin::FeatureSet
AEIaccPlugin::process(const float *const *inputBuffers, Vamp::RealTime timestamp)
{
    FeatureSet fs;
    (void)timestamp;
    
    // Process each channel independently
    for (size_t ch = 0; ch < m_channels; ++ch) {
        processBlock(ch, inputBuffers[ch]);
    }
    
    // The counts that need the global maximum are made in getRemainingFeatures()
    return fs;
}

//Bnakers' rounding to nearest integer
static int round_r(double x) {
    double rounded = std::floor(x + 0.5);
    // If exactly at .5, round to nearest even
    if (x - std::floor(x) == 0.5) {
        int floored = static_cast<int>(std::floor(x));
        if (floored % 2 != 0) {
            rounded = std::floor(x);
        }
    }
    return static_cast<int>(rounded);
}

// Set up the bands (same for all channels); done once, when the results are first needed
void AEIaccPlugin::initBands()
{
    if (m_bandsInitialized) return;

    float maxFreqHz = m_maxFreq;
    if (maxFreqHz <= 0) maxFreqHz = m_inputSampleRate / 2.0f;
    
    int numBands = static_cast<int>(std::ceil(maxFreqHz / m_freqStep));
    
    // Initialize per-channel band counters
    for (size_t ch = 0; ch < m_channels; ++ch) {
        m_bandAboveThresholdCount_ch[ch].assign(numBands, 0);
        m_bandTotalCount_ch[ch].assign(numBands, 0);
    }
    
    m_bandStartBins.assign(numBands, -1);
    m_bandEndBins.assign(numBands, -1);
    
    // Calculate band boundaries matching soundecology logic
    for (int i = 0; i < numBands; ++i) {
        float minFreq = i * m_freqStep;
        float maxFreq = (i + 1) * m_freqStep;
        
        int miny = round_r(minFreq / 10.0);
        int maxy = round_r(maxFreq / 10.0);
        
        int startBin = (miny <= 1) ? 0 : miny - 1;
        int endBin = maxy;
        
        if (startBin < 0) startBin = 0;
        if (endBin > static_cast<int>(m_numBins)) endBin = static_cast<int>(m_numBins);
        
        if (startBin < endBin) {
            m_bandStartBins[i] = startBin;
            m_bandEndBins[i] = endBin;
        }
    }
    m_bandsInitialized = true;
}

// Transform one block of one channel and keep its scaled magnitudes. This is what both passes over the frames
// used to compute (the first only for the global maximum, the second again for the threshold counts); the
// values and the order in which the maximum is updated are the same, the FFT is done once.
void AEIaccPlugin::processBlock(size_t channel, const float* input)
{
    const size_t blockSize = m_blockSize;
    double* frameData = m_fft->buffer();
    for (size_t i = 0; i < blockSize; ++i) {
        frameData[i] = input[i] * m_window[i];
    }
    m_fft->forward();

    double* magnitudes = m_magnitude_ch[channel].appendRow();
    float& globalMax = m_globalMax_ch[channel];
    const size_t windowLength = m_blockSize;
    m_fft->forBins(0, m_numBins, [magnitudes, &globalMax, windowLength](size_t b, double real, double imag) {
        double magnitude = std::sqrt(real * real + imag * imag);
        // Apply R's scaling: divide by window length, multiply by 2
        magnitude = magnitude / windowLength * 2.0;
        if (static_cast<float>(magnitude) > globalMax) globalMax = static_cast<float>(magnitude);
        magnitudes[b] = magnitude;
    });
    
    m_frameCount_ch[channel]++;
}

Vamp::Plugin::FeatureSet
AEIaccPlugin::getRemainingFeatures()
{
    FeatureSet fs;
    
    // The global max for each channel is known: every block has been transformed
    // Calculate AEI for each channel
    for (size_t ch = 0; ch < m_channels; ++ch) {
        if (m_frameCount_ch[ch] == 0) {
            Feature f;
            f.hasTimestamp = true;
            f.timestamp = Vamp::RealTime::frame2RealTime(0, static_cast<unsigned int>(m_inputSampleRate));
            f.values.push_back(0.0f);
            fs[static_cast<int>(ch)].push_back(f);
            continue;
        }

        initBands();

        // Calculate global max dB and threshold for this channel
        double globalMaxDB = -200.0;
        if (m_globalMax_ch[ch] > 1e-10) {
            globalMaxDB = 20.0 * std::log10(m_globalMax_ch[ch]);
        }
        
        double thresholdDB = globalMaxDB + m_dbThreshold;
        
        // Count values above threshold for each band. A cell is above the threshold when
        // 20 * log10(magnitude) > thresholdDB (-200 dB for magnitudes up to 1e-10); DbAbove gives that answer and
        // calls log10 only for magnitudes close to the threshold (see DbAbove.h).
        size_t numFrames = m_magnitude_ch[ch].rows();
        size_t numBands = m_bandAboveThresholdCount_ch[ch].size();
        const DbAbove isAbove(thresholdDB);

        for (size_t frame = 0; frame < numFrames; ++frame) {
            const double* magnitudes = m_magnitude_ch[ch].row(frame);

            for (size_t band = 0; band < numBands; ++band) {
                int start = m_bandStartBins[band];
                int end = m_bandEndBins[band];

                if (start != -1 && end != -1) {
                    long above = 0;
                    for (int b = start; b < end; ++b) {
                        if (isAbove(magnitudes[b])) {
                            ++above;
                        }
                    }

                    // Count total cells and cells above threshold
                    m_bandTotalCount_ch[ch][band] += end - start;
                    m_bandAboveThresholdCount_ch[ch][band] += above;
                }
            }
        }

        m_magnitude_ch[ch].clear();
        
        // Calculate band proportions (cells above threshold / total cells)
        std::vector<double> bandCounts(numBands, 0.0);
        
        for (size_t i = 0; i < numBands; ++i) {
            if (m_bandTotalCount_ch[ch][i] > 0) {
                bandCounts[i] = static_cast<double>(m_bandAboveThresholdCount_ch[ch][i]) / 
                               static_cast<double>(m_bandTotalCount_ch[ch][i]);
            } else {
                bandCounts[i] = 0.0;
            }
        }
        
        // Calculate Gini Index
        std::sort(bandCounts.begin(), bandCounts.end());
        
        double sumCounts = 0.0;
        double weightedSum = 0.0;
        double n = static_cast<double>(numBands);
        
        for (size_t i = 0; i < numBands; ++i) {
            sumCounts += bandCounts[i];
            weightedSum += (i + 1) * bandCounts[i];
        }
        
        double aei = 0.0;
        if (sumCounts > 0) {
            aei = (2.0 * weightedSum) / (n * sumCounts) - (n + 1.0) / n;
        }
        
        Feature f;
        f.hasTimestamp = true;
        f.timestamp = Vamp::RealTime::frame2RealTime(0, static_cast<unsigned int>(m_inputSampleRate));
        f.values.push_back(static_cast<float>(aei));
        
        fs[static_cast<int>(ch)].push_back(f);
    }
    
    return fs;
}
