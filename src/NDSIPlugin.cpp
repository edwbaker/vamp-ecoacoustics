/*
    Ecoacoustic Vamp Plugins
    
    High-performance implementations of acoustic indices for 
    bioacoustics and soundscape ecology analysis.
    
    (C) Ed Baker 2025. Licensed under GPL (>=2).
*/

// NDSIPlugin.cpp - Normalized Difference Soundscape Index
// Matches soundecology::ndsi behavior

#include "NDSIPlugin.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <cstdio>
#include <limits>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using std::string;
using std::vector;

NDSIPlugin::NDSIPlugin(float inputSampleRate) :
    Plugin(inputSampleRate),
    m_channels(0),
    m_blockSize(0),
    m_stepSize(0),
    m_fftSize(1024),
    m_numBins(512),
    m_secondCount(0),
    m_anthroMin(1000.0f),
    m_anthroMax(2000.0f),
    m_bioMin(2000.0f),
    m_bioMax(11000.0f),
    m_fft(nullptr),
    m_windowSumSq(0.0)
{
}

NDSIPlugin::~NDSIPlugin()
{
    if (m_fft) {
        delete m_fft;
        m_fft = nullptr;
    }
}

string
NDSIPlugin::getIdentifier() const
{
    return "ndsi";
}

string
NDSIPlugin::getName() const
{
    return "Normalized Difference Soundscape Index (NDSI)";
}

string
NDSIPlugin::getDescription() const
{
    return "Calculates the Normalized Difference Soundscape Index (NDSI) matching soundecology::ndsi.";
}

string
NDSIPlugin::getMaker() const
{
    return "Ecoacoustic-Vamp-Plugins";
}

int
NDSIPlugin::getPluginVersion() const
{
    return 2;
}

string
NDSIPlugin::getCopyright() const
{
    return "(C) Ed Baker 2025. Licensed under GPL (>=2).";
}

Vamp::Plugin::InputDomain
NDSIPlugin::getInputDomain() const
{
    return TimeDomain;
}

size_t
NDSIPlugin::getMinChannelCount() const
{
    return 1;
}

size_t
NDSIPlugin::getMaxChannelCount() const
{
    return 2;
}

NDSIPlugin::ParameterList NDSIPlugin::getParameterDescriptors() const
{
    ParameterList list;

    ParameterDescriptor d;
    d.identifier = "anthroMin";
    d.name = "Anthropophony Min Freq";
    d.description = "Minimum frequency for anthropophony (Hz)";
    d.unit = "Hz";
    d.minValue = 0;
    d.maxValue = m_inputSampleRate / 2;
    d.defaultValue = 1000;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "anthroMax";
    d.name = "Anthropophony Max Freq";
    d.description = "Maximum frequency for anthropophony (Hz)";
    d.defaultValue = 2000;
    list.push_back(d);

    d.identifier = "bioMin";
    d.name = "Biophony Min Freq";
    d.description = "Minimum frequency for biophony (Hz)";
    d.defaultValue = 2000;
    list.push_back(d);

    d.identifier = "bioMax";
    d.name = "Biophony Max Freq";
    d.description = "Maximum frequency for biophony (Hz)";
    d.defaultValue = 11000;
    list.push_back(d);

    return list;
}

float NDSIPlugin::getParameter(string identifier) const
{
    if (identifier == "anthroMin") return m_anthroMin;
    if (identifier == "anthroMax") return m_anthroMax;
    if (identifier == "bioMin") return m_bioMin;
    if (identifier == "bioMax") return m_bioMax;
    return 0;
}

void NDSIPlugin::setParameter(string identifier, float value)
{
    if (identifier == "anthroMin") m_anthroMin = value;
    if (identifier == "anthroMax") m_anthroMax = value;
    if (identifier == "bioMin") m_bioMin = value;
    if (identifier == "bioMax") m_bioMax = value;
}

NDSIPlugin::ProgramList
NDSIPlugin::getPrograms() const
{
    return ProgramList();
}

string
NDSIPlugin::getCurrentProgram() const
{
    return "";
}

void
NDSIPlugin::selectProgram(string program)
{
    (void)program;
}

NDSIPlugin::OutputList NDSIPlugin::getOutputDescriptors() const
{
    OutputList list;

    OutputDescriptor d;
    d.identifier = "ndsi";
    d.name = "NDSI";
    d.description = "Normalized Difference Soundscape Index";
    d.unit = "";
    d.hasFixedBinCount = true;
    d.binCount = 1;
    d.hasKnownExtents = true;
    d.minValue = -1;
    d.maxValue = 1;
    d.isQuantized = false;
    d.sampleType = OutputDescriptor::VariableSampleRate;
    d.hasDuration = false;
    list.push_back(d);

    OutputDescriptor d2;
    d2.identifier = "biophony";
    d2.name = "Biophony";
    d2.description = "Biophony Energy (normalized)";
    d2.sampleType = OutputDescriptor::VariableSampleRate;
    d2.hasFixedBinCount = true;
    d2.binCount = 1;
    list.push_back(d2);

    OutputDescriptor d3;
    d3.identifier = "anthropophony";
    d3.name = "Anthropophony";
    d3.description = "Anthropophony Energy (normalized)";
    d3.sampleType = OutputDescriptor::VariableSampleRate;
    d3.hasFixedBinCount = true;
    d3.binCount = 1;
    list.push_back(d3);

    // Right channel outputs (for stereo files)
    OutputDescriptor d4;
    d4.identifier = "ndsi_right";
    d4.name = "NDSI (Right)";
    d4.description = "Normalized Difference Soundscape Index - Right Channel";
    d4.unit = "";
    d4.hasFixedBinCount = true;
    d4.binCount = 1;
    d4.hasKnownExtents = true;
    d4.minValue = -1;
    d4.maxValue = 1;
    d4.isQuantized = false;
    d4.sampleType = OutputDescriptor::VariableSampleRate;
    d4.hasDuration = false;
    list.push_back(d4);

    OutputDescriptor d5;
    d5.identifier = "biophony_right";
    d5.name = "Biophony (Right)";
    d5.description = "Biophony Energy (normalized) - Right Channel";
    d5.sampleType = OutputDescriptor::VariableSampleRate;
    d5.hasFixedBinCount = true;
    d5.binCount = 1;
    list.push_back(d5);

    OutputDescriptor d6;
    d6.identifier = "anthropophony_right";
    d6.name = "Anthropophony (Right)";
    d6.description = "Anthropophony Energy (normalized) - Right Channel";
    d6.sampleType = OutputDescriptor::VariableSampleRate;
    d6.hasFixedBinCount = true;
    d6.binCount = 1;
    list.push_back(d6);

    return list;
}

bool NDSIPlugin::initialise(size_t channels, size_t stepSize, size_t blockSize)
{
    if (channels < getMinChannelCount() || channels > getMaxChannelCount()) return false;

    m_channels = channels;
    m_blockSize = blockSize;
    m_stepSize = stepSize;
    
    m_fftSize = 1024;
    m_numBins = m_fftSize / 2;
    
    if (m_fft) delete m_fft;
    m_fft = new Vamp::FFTReal(static_cast<unsigned int>(m_fftSize));
    m_fftOut.resize(m_fftSize + 2);
    
    // Hamming window
    m_window.resize(m_fftSize);
    m_windowSumSq = 0.0;
    for (size_t i = 0; i < m_fftSize; ++i) {
        double w = 0.54 - 0.46 * std::cos(2.0 * M_PI * i / (m_fftSize - 1));
        m_window[i] = static_cast<float>(w);
        m_windowSumSq += w * w;
    }
    
    // Pre-allocate work buffers for processOneSecond
    m_timeBuf.resize(m_fftSize);
    m_detrendBuf.resize(m_fftSize);
    m_windowPSD.resize(m_numBins);
    
    // Pre-compute detrending constants (sumT = 1+2+...+N, sumTT = 1²+2²+...+N²)
    double nD = static_cast<double>(m_fftSize);
    m_detrendSumT = nD * (nD + 1.0) / 2.0;  // Sum of 1..N
    m_detrendSumTT = nD * (nD + 1.0) * (2.0 * nD + 1.0) / 6.0;  // Sum of 1²..N²
    m_detrendDenom = nD * m_detrendSumTT - m_detrendSumT * m_detrendSumT;
    m_detrendNfftInv = 1.0 / nD;
    
    // Pre-compute PSD normalization factors
    m_psdNorm = 2.0 / (m_windowSumSq * m_fftSize);
    m_psdNyquistNorm = 1.0 / (m_windowSumSq * m_fftSize);

    // Initialize per-channel data structures
    m_sampleBuffer.resize(m_channels);
    m_accumulatedPSD.resize(m_channels);
    m_secondCount.resize(m_channels, 0);
    
    for (size_t ch = 0; ch < m_channels; ++ch) {
        m_sampleBuffer[ch].clear();
        m_accumulatedPSD[ch].assign(m_numBins, 0.0);
    }
    
    // Initialize padding detection
    m_totalSamplesReceived = 0;
    m_lastBlockValidSamples = 0;

    return true;
}

void NDSIPlugin::reset()
{
    for (size_t ch = 0; ch < m_channels; ++ch) {
        m_secondCount[ch] = 0;
        std::fill(m_accumulatedPSD[ch].begin(), m_accumulatedPSD[ch].end(), 0.0);
        m_sampleBuffer[ch].clear();
    }
    m_totalSamplesReceived = 0;
    m_lastBlockValidSamples = 0;
}

// Process one second of audio using Welch's method
void NDSIPlugin::processOneSecond(const vector<double>& samples, size_t channel)
{
    const size_t n = samples.size();
    const size_t nfft = m_fftSize;
    const size_t noverlap = nfft / 2;
    
    // Calculate number of windows
    const size_t pwelchStep = nfft - noverlap + 1;
    
    size_t start = 0;
    size_t nWindows = 0;
    
    // Zero out reusable work buffers
    std::fill(m_windowPSD.begin(), m_windowPSD.end(), 0.0);
    
    const float* win = m_window.data();
    const double* sampleData = samples.data();
    
    while (start + nfft <= n) {
        const double* src = sampleData + start;
        for (size_t i = 0; i < nfft; ++i) {
            m_timeBuf[i] = src[i] * win[i];
        }
        
        // Detrend using pre-computed constants
        double sumY = 0.0, sumTY = 0.0;
        for (size_t i = 0; i < nfft; ++i) {
            double t = static_cast<double>(i + 1);
            double y = m_timeBuf[i];
            sumY += y;
            sumTY += t * y;
        }
        double b = (m_detrendDenom != 0.0) ? (static_cast<double>(nfft) * sumTY - m_detrendSumT * sumY) / m_detrendDenom : 0.0;
        double a = (sumY - b * m_detrendSumT) * m_detrendNfftInv;
        for (size_t i = 0; i < nfft; ++i) {
            double t = static_cast<double>(i + 1);
            m_detrendBuf[i] = m_timeBuf[i] - (a + b * t);
        }
        
        // Zero pad fftOut remainder
        for (size_t i = nfft; i < m_fftOut.size(); ++i) m_fftOut[i] = 0.0;

        // Execute FFT on detrended windowed data
        m_fft->forward(m_detrendBuf.data(), m_fftOut.data());
        
        // Compute one-sided PSD for positive frequencies (skip DC)
        const double* fftOut = m_fftOut.data();
        for (size_t i = 1; i < m_numBins; ++i) {
            double real = fftOut[2 * i];
            double imag = fftOut[2 * i + 1];
            m_windowPSD[i - 1] += (real * real + imag * imag) * m_psdNorm;
        }
        // Nyquist bin
        {
            double real = fftOut[2 * m_numBins];
            double imag = fftOut[2 * m_numBins + 1];
            m_windowPSD[m_numBins - 1] += (real * real + imag * imag) * m_psdNyquistNorm;
        }
        
        nWindows++;
        start += pwelchStep;
    }
    
    if (nWindows == 0) return;

    // Average over windows and add to accumulated PSD for this channel
    const double invWindows = 1.0 / static_cast<double>(nWindows);
    double* accPSD = m_accumulatedPSD[channel].data();
    for (size_t i = 0; i < m_numBins; ++i) {
        accPSD[i] += m_windowPSD[i] * invWindows;
    }
    
    m_secondCount[channel]++;
}

Vamp::Plugin::FeatureSet NDSIPlugin::process(const float *const *inputBuffers, Vamp::RealTime timestamp)
{
    FeatureSet fs;
    (void)timestamp;
    
    // Detect trailing zeros (padding) in this block
    size_t validSamples = m_blockSize;
    for (size_t i = m_blockSize; i > 0; --i) {
        if (inputBuffers[0][i-1] != 0.0f) {
            validSamples = i;
            break;
        }
        if (i == 1) {
            validSamples = 0;
        }
    }
    
    // Track samples for padding detection
    m_totalSamplesReceived += m_blockSize;
    m_lastBlockValidSamples = validSamples;
    
    // Buffer incoming samples for each channel separately
    for (size_t ch = 0; ch < m_channels; ++ch) {
        size_t oldSize = m_sampleBuffer[ch].size();
        m_sampleBuffer[ch].resize(oldSize + m_blockSize);
        double* dst = m_sampleBuffer[ch].data() + oldSize;
        const float* src = inputBuffers[ch];
        for (size_t i = 0; i < m_blockSize; ++i) {
            dst[i] = static_cast<double>(src[i]);
        }
    }
    
    // Process complete 1-second chunks, but KEEP ONE SECOND buffered
    // This ensures we don't process a second that might contain padding
    size_t samplesPerSecond = static_cast<size_t>(m_inputSampleRate);
    
    while (m_sampleBuffer[0].size() >= 2 * samplesPerSecond) {
        for (size_t ch = 0; ch < m_channels; ++ch) {
            vector<double> oneSecond(m_sampleBuffer[ch].begin(), m_sampleBuffer[ch].begin() + samplesPerSecond);
            processOneSecond(oneSecond, ch);
            m_sampleBuffer[ch].erase(m_sampleBuffer[ch].begin(), m_sampleBuffer[ch].begin() + samplesPerSecond);
        }
    }
    
    return fs;
}

size_t
NDSIPlugin::getPreferredBlockSize() const
{
    return 1024;
}

size_t
NDSIPlugin::getPreferredStepSize() const
{
    return 1024;
}

// Calculate NDSI for a specific channel from accumulated PSD
void NDSIPlugin::calculateNDSI(size_t channel, double& ndsi, double& biophony, double& anthrophony)
{
    ndsi = std::numeric_limits<double>::quiet_NaN();
    biophony = std::numeric_limits<double>::quiet_NaN();
    anthrophony = std::numeric_limits<double>::quiet_NaN();
    
    if (m_secondCount[channel] == 0) return;
    
    // Average the accumulated PSD for this channel
    vector<double> meanPSD = m_accumulatedPSD[channel];
    for (size_t i = 0; i < meanPSD.size(); ++i) {
        meanPSD[i] /= static_cast<double>(m_secondCount[channel]);
    }

    const double nyquist = m_inputSampleRate / 2.0;
    const int specRows = static_cast<int>(meanPSD.size());
    if (specRows <= 0 || nyquist <= 0.0) {
        return;
    }

    const double freqPerRow = static_cast<double>(specRows) / nyquist;
    double hzInterval = static_cast<double>(m_anthroMax) - static_cast<double>(m_anthroMin);
    if (hzInterval <= 0.0) {
        hzInterval = 1.0;
    }

    auto clampRow = [specRows](double row) -> int {
        int idx = static_cast<int>(std::round(row));
        if (idx < 1) idx = 1;
        if (idx > specRows) idx = specRows;
        return idx;
    };

    // trapzRange: compute trapezoidal integral over 1-based indices [startRow1, endRow1]
    auto trapzRange = [&meanPSD, specRows](int startRow1, int endRow1) -> double {
        if (endRow1 < startRow1) {
            return 0.0;
        }

        int startIdx = startRow1 - 1;
        int endIdx = endRow1 - 1;  // This is now the LAST element index (inclusive)
        if (startIdx < 0) startIdx = 0;
        if (endIdx >= specRows) endIdx = specRows - 1;
        if (endIdx < startIdx) {
            return 0.0;
        }

        double sum = 0.0;
        for (int i = startIdx; i < endIdx; ++i) {
            sum += (meanPSD[static_cast<size_t>(i)] + meanPSD[static_cast<size_t>(i + 1)]) * 0.5;
        }
        return sum;
    };

    const double anthroRangeHz = static_cast<double>(m_anthroMax) - static_cast<double>(m_anthroMin);
    int anthroBinsCount = static_cast<int>(std::round(anthroRangeHz / hzInterval));
    if (anthroBinsCount < 1) anthroBinsCount = 1;
    const int bioBinsCount = std::max(1, static_cast<int>(std::round((static_cast<double>(m_bioMax) - static_cast<double>(m_bioMin)) / hzInterval)));

    const int anthroMinRow = clampRow(static_cast<double>(m_anthroMin) * freqPerRow);
    const int anthroMaxRow = clampRow(static_cast<double>(m_anthroMax) * freqPerRow);

    std::vector<double> anthroBins(static_cast<size_t>(anthroBinsCount), 0.0);
    for (int i = 0; i < anthroBinsCount; ++i) {
        anthroBins[static_cast<size_t>(i)] = trapzRange(anthroMinRow, anthroMaxRow);
    }

    const double bioRangeHz = static_cast<double>(m_bioMax) - static_cast<double>(m_bioMin);
    const double bioStepRange = (bioBinsCount > 0) ? freqPerRow * (bioRangeHz / static_cast<double>(bioBinsCount)) : 0.0;
    double bioMinRowDouble = static_cast<double>(std::round(static_cast<double>(m_bioMin) * freqPerRow));
    double bioMaxRowDouble = bioMinRowDouble + bioStepRange;
    std::vector<double> bioBins(static_cast<size_t>(bioBinsCount), 0.0);
    
    for (int i = 0; i < bioBinsCount; ++i) {
        int startRow = static_cast<int>(std::floor(bioMinRowDouble));
        int endRow = startRow + static_cast<int>(std::floor(bioMaxRowDouble - bioMinRowDouble));
        
        if (startRow < 1) startRow = 1;
        if (endRow >= specRows) endRow = specRows - 1;
        if (endRow < startRow) endRow = startRow;
        
        bioBins[static_cast<size_t>(i)] = trapzRange(startRow, endRow);
        bioMinRowDouble += bioStepRange;
        bioMaxRowDouble += bioStepRange;
    }

    std::vector<double> freqBins;
    freqBins.reserve(anthroBins.size() + bioBins.size());
    freqBins.insert(freqBins.end(), anthroBins.begin(), anthroBins.end());
    freqBins.insert(freqBins.end(), bioBins.begin(), bioBins.end());

    double frob = std::sqrt(std::inner_product(freqBins.begin(), freqBins.end(), freqBins.begin(), 0.0));
    if (frob > 0.0) {
        for (double &value : freqBins) {
            value /= frob;
        }
    } else if (!freqBins.empty()) {
        for (double &value : freqBins) {
            value = std::numeric_limits<double>::quiet_NaN();
        }
    }

    anthrophony = freqBins.empty() ? 0.0 : freqBins.front();
    biophony = 0.0;
    if (freqBins.size() > 1) {
        for (size_t i = 1; i < freqBins.size(); ++i) {
            biophony += freqBins[i];
        }
    }

    double denom = biophony + anthrophony;
    if (denom != 0.0) {
        ndsi = (biophony - anthrophony) / denom;
    } else {
        ndsi = std::numeric_limits<double>::quiet_NaN();
    }
}

NDSIPlugin::FeatureSet NDSIPlugin::getRemainingFeatures()
{
    FeatureSet fs;
    
    // Calculate true sample count (excluding padding)
    size_t paddingSamples = m_blockSize - m_lastBlockValidSamples;
    size_t trueSampleCount = m_totalSamplesReceived - paddingSamples;
    
    // Calculate how many complete seconds we should have
    double trueDuration = static_cast<double>(trueSampleCount) / m_inputSampleRate;
    size_t expectedSeconds = static_cast<size_t>(std::floor(trueDuration));
    
    // Process any remaining complete seconds from the buffer
    size_t samplesPerSecond = static_cast<size_t>(m_inputSampleRate);
    
    // Process remaining complete seconds for each channel
    while (m_secondCount[0] < expectedSeconds && m_sampleBuffer[0].size() >= samplesPerSecond) {
        for (size_t ch = 0; ch < m_channels; ++ch) {
            vector<double> oneSecond(m_sampleBuffer[ch].begin(), m_sampleBuffer[ch].begin() + samplesPerSecond);
            processOneSecond(oneSecond, ch);
            m_sampleBuffer[ch].erase(m_sampleBuffer[ch].begin(), m_sampleBuffer[ch].begin() + samplesPerSecond);
        }
    }

    // Calculate NDSI for left channel (or mono)
    double ndsiLeft, bioLeft, anthroLeft;
    calculateNDSI(0, ndsiLeft, bioLeft, anthroLeft);
    
    if (m_secondCount[0] == 0) return fs;
    
    // Output left/mono channel results
    Feature f;
    f.hasTimestamp = true;
    f.timestamp = Vamp::RealTime::zeroTime;
    f.values.push_back(static_cast<float>(ndsiLeft));
    fs[0].push_back(f);
    
    Feature fB;
    fB.hasTimestamp = true;
    fB.timestamp = Vamp::RealTime::zeroTime;
    fB.values.push_back(static_cast<float>(bioLeft));
    fs[1].push_back(fB);
    
    Feature fA;
    fA.hasTimestamp = true;
    fA.timestamp = Vamp::RealTime::zeroTime;
    fA.values.push_back(static_cast<float>(anthroLeft));
    fs[2].push_back(fA);
    
    // Calculate and output right channel results if stereo
    if (m_channels == 2) {
        double ndsiRight, bioRight, anthroRight;
        calculateNDSI(1, ndsiRight, bioRight, anthroRight);
        
        Feature fR;
        fR.hasTimestamp = true;
        fR.timestamp = Vamp::RealTime::zeroTime;
        fR.values.push_back(static_cast<float>(ndsiRight));
        fs[3].push_back(fR);
        
        Feature fBR;
        fBR.hasTimestamp = true;
        fBR.timestamp = Vamp::RealTime::zeroTime;
        fBR.values.push_back(static_cast<float>(bioRight));
        fs[4].push_back(fBR);
        
        Feature fAR;
        fAR.hasTimestamp = true;
        fAR.timestamp = Vamp::RealTime::zeroTime;
        fAR.values.push_back(static_cast<float>(anthroRight));
        fs[5].push_back(fAR);
    }
    
    return fs;
}
