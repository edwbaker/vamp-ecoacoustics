// ADIaccPlugin.cpp - Acoustic Diversity Index (Accumulated)
// Matches soundecology::acoustic_diversity behavior

#include "ADIaccPlugin.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <map>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

ADIaccPlugin::ADIaccPlugin(float inputSampleRate) :
    EcoacousticSpectralPlugin(inputSampleRate),
    m_dbThreshold(-50.0f),
    m_freqStep(1000.0f)
{
    m_maxFreq = 10000.0f; // Default max freq for ADI
}

ADIaccPlugin::~ADIaccPlugin()
{
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
    return "ReVAMP";
}

int
ADIaccPlugin::getPluginVersion() const
{
    return 1;
}

string
ADIaccPlugin::getCopyright() const
{
    return "MIT License";
}

Vamp::Plugin::ParameterList
ADIaccPlugin::getParameterDescriptors() const
{
    ParameterList list;

    ParameterDescriptor d;
    d.identifier = "maxfreq";
    d.name = "Maximum Frequency";
    d.description = "Maximum frequency (kHz) to include in calculation";
    d.unit = "kHz";
    d.minValue = 0;
    d.maxValue = m_inputSampleRate / 2000.0f;
    d.defaultValue = 10.0f;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "db_threshold";
    d.name = "dB Threshold";
    d.description = "Threshold in dBFS";
    d.unit = "dB";
    d.minValue = -120.0f;
    d.maxValue = 0.0f;
    d.defaultValue = -50.0f;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "freq_step";
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
    if (identifier == "maxfreq") return m_maxFreq; // Stored in base class but used here
    if (identifier == "db_threshold") return m_dbThreshold;
    if (identifier == "freq_step") return m_freqStep;
    return 0;
}

void
ADIaccPlugin::setParameter(string identifier, float value)
{
    if (identifier == "maxfreq") {
        m_maxFreq = value;
    } else if (identifier == "db_threshold") {
        m_dbThreshold = value;
    } else if (identifier == "freq_step") {
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
}

Vamp::Plugin::OutputList
ADIaccPlugin::getOutputDescriptors() const
{
    OutputList list;

    OutputDescriptor d;
    d.identifier = "adi";
    d.name = "ADI";
    d.description = "Acoustic Diversity Index";
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

size_t
ADIaccPlugin::getPreferredBlockSize() const
{
    return 4096;
}

size_t
ADIaccPlugin::getPreferredStepSize() const
{
    return 4096;
}

bool
ADIaccPlugin::initialise(size_t channels, size_t stepSize, size_t blockSize)
{
    if (!EcoacousticSpectralPlugin::initialise(channels, stepSize, blockSize)) return false;

    m_bandsInitialized = false;
    m_bandHistograms.clear();
    m_bandStartBins.clear();
    m_bandEndBins.clear();

    return true;
}

void ADIaccPlugin::processBatch(size_t numFrames)
{
    if (numFrames == 0) return;

    // Initialize bands if needed
    if (!m_bandsInitialized) {
        float maxFreqHz = m_maxFreq;
        if (maxFreqHz <= 0) maxFreqHz = m_inputSampleRate / 2.0f;
        
        int numBands = static_cast<int>(std::ceil(maxFreqHz / m_freqStep));
        
        // 3000 bins for -200.0 to +100.0 dB with 0.1 step
        m_bandHistograms.assign(numBands, std::vector<int>(3000, 0));
        
        m_bandStartBins.assign(numBands, -1);
        m_bandEndBins.assign(numBands, -1);
        
        float binResolution = (m_inputSampleRate / 2.0f) / (m_blockSize / 2);
        
        for (size_t b = 0; b < m_numBins; ++b) {
            float freq = (b + 1) * binResolution;
            if (freq <= maxFreqHz) {
                int bandIndex = static_cast<int>(std::floor(freq / m_freqStep));
                if (bandIndex >= 0 && bandIndex < numBands) {
                    if (m_bandStartBins[bandIndex] == -1) {
                        m_bandStartBins[bandIndex] = b;
                    }
                    m_bandEndBins[bandIndex] = b + 1;
                }
            }
        }
        m_bandsInitialized = true;
    }

    size_t blockSize = m_blockSize;
    size_t numBands = m_bandHistograms.size();

    for (size_t frame = 0; frame < numFrames; ++frame) {
        double* frameData = &m_inputBuffer[frame * blockSize];
        m_fft->forward(frameData, m_fftOut.data());
        
        // Optimized loop: Iterate bands
        for (size_t band = 0; band < numBands; ++band) {
            int start = m_bandStartBins[band];
            int end = m_bandEndBins[band];
            
            if (start != -1 && end != -1) {
                for (int b = start; b < end; ++b) {
                    // FFT bin index is b+1
                    double real = m_fftOut[2 * (b + 1)];
                    double imag = m_fftOut[2 * (b + 1) + 1];
                    double magnitude = std::sqrt(real * real + imag * imag);
                    
                    if (magnitude > m_globalMax) m_globalMax = magnitude;
                    
                    double db = -200.0;
                    if (magnitude > 1e-10) {
                        db = 20.0 * std::log10(magnitude);
                    }
                    
                    // Map to bin: -200 to +100 dB relative to 1.0
                    int bin = static_cast<int>((db + 200.0) * 10.0);
                    if (bin < 0) bin = 0;
                    if (bin >= 3000) bin = 2999;
                    
                    m_bandHistograms[band][bin]++;
                }
            }
        }
        
        m_frameCount++;
    }
}

Vamp::Plugin::FeatureSet
ADIaccPlugin::getRemainingFeatures()
{
    FeatureSet fs;
    
    if (!m_inputBuffer.empty()) {
        size_t remainingFrames = m_inputBuffer.size() / m_blockSize;
        if (remainingFrames > 0) {
            processBatch(remainingFrames);
        }
        m_inputBuffer.clear();
    }
    
    if (m_frameCount == 0) return fs;

    // Calculate Threshold relative to Global Max
    double globalMaxDB = -200.0;
    if (m_globalMax > 1e-10) {
        globalMaxDB = 20.0 * std::log10(m_globalMax);
    }
    
    double thresholdDB = globalMaxDB + m_dbThreshold;
    
    // Convert threshold to histogram bin
    int thresholdBin = static_cast<int>((thresholdDB + 200.0) * 10.0);
    
    // Calculate band proportions
    size_t numBands = m_bandHistograms.size();
    std::vector<double> bandProportions(numBands, 0.0);
    double sumProportions = 0.0;
    
    for (size_t i = 0; i < numBands; ++i) {
        long count = 0;
        long totalCells = 0;
        
        for (size_t bin = 0; bin < m_bandHistograms[i].size(); ++bin) {
            totalCells += m_bandHistograms[i][bin];
            if (static_cast<int>(bin) > thresholdBin) {
                count += m_bandHistograms[i][bin];
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
                adi -= p * std::log(p);
            }
        }
    }
    
    Feature f;
    f.hasTimestamp = true;
    f.timestamp = Vamp::RealTime::frame2RealTime(0, m_inputSampleRate);
    f.values.push_back(static_cast<float>(adi));
    
    fs[0].push_back(f);
    
    return fs;
}
