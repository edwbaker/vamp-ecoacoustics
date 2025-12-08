#include "BIPlugin.h"
#include <cmath>
#include <iostream>

BIPlugin::BIPlugin(float inputSampleRate) :
    EcoacousticSpectralPlugin(inputSampleRate),
    m_minFreq(2000.0f),
    m_maxFreq(15000.0f)
{
}

BIPlugin::~BIPlugin()
{
}

string BIPlugin::getIdentifier() const
{
    return "bi";
}

string BIPlugin::getName() const
{
    return "Bioacoustic Index (Maad)";
}

string BIPlugin::getDescription() const
{
    return "Bioacoustic Index calculated from the mean spectrum in dB (Scikit-Maad style)";
}

string BIPlugin::getMaker() const
{
    return "Your Name";
}

int BIPlugin::getPluginVersion() const
{
    return 1;
}

string BIPlugin::getCopyright() const
{
    return "MIT";
}

BIPlugin::ParameterList BIPlugin::getParameterDescriptors() const
{
    ParameterList list;

    ParameterDescriptor d;
    d.identifier = "min_freq";
    d.name = "Minimum Frequency";
    d.description = "Lower bound of the frequency band (Hz)";
    d.unit = "Hz";
    d.minValue = 0;
    d.maxValue = m_inputSampleRate / 2;
    d.defaultValue = 2000.0f;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "max_freq";
    d.name = "Maximum Frequency";
    d.description = "Upper bound of the frequency band (Hz)";
    d.unit = "Hz";
    d.minValue = 0;
    d.maxValue = m_inputSampleRate / 2;
    d.defaultValue = 15000.0f;
    d.isQuantized = false;
    list.push_back(d);

    return list;
}

float BIPlugin::getParameter(string identifier) const
{
    if (identifier == "min_freq") return m_minFreq;
    if (identifier == "max_freq") return m_maxFreq;
    return 0;
}

void BIPlugin::setParameter(string identifier, float value)
{
    if (identifier == "min_freq") m_minFreq = value;
    if (identifier == "max_freq") m_maxFreq = value;
}

BIPlugin::ProgramList BIPlugin::getPrograms() const
{
    ProgramList list;
    return list;
}

string BIPlugin::getCurrentProgram() const
{
    return "";
}

void BIPlugin::selectProgram(string name)
{
}

BIPlugin::OutputList BIPlugin::getOutputDescriptors() const
{
    OutputList list;

    OutputDescriptor d;
    d.identifier = "bi";
    d.name = "Bioacoustic Index";
    d.description = "The Bioacoustic Index (Maad)";
    d.unit = "";
    d.hasFixedBinCount = true;
    d.binCount = 1;
    d.hasKnownExtents = false;
    d.isQuantized = false;
    d.sampleType = OutputDescriptor::OneSamplePerStep; // Actually OneSamplePerFile, but we return at end
    list.push_back(d);

    return list;
}

bool BIPlugin::initialise(size_t channels, size_t stepSize, size_t blockSize)
{
    if (!EcoacousticSpectralPlugin::initialise(channels, stepSize, blockSize)) return false;
    
    m_accumulatedSpectrum.assign(m_numBins, 0.0);
    
    return true;
}

size_t BIPlugin::getPreferredBlockSize() const
{
    return 512;
}

size_t BIPlugin::getPreferredStepSize() const
{
    return 512; // No overlap
}

void BIPlugin::processBatch(size_t numFrames)
{
    if (numFrames == 0) return;

    size_t blockSize = m_blockSize;
    
    for (size_t frame = 0; frame < numFrames; ++frame) {
        // Get pointer to current frame in input buffer
        double* frameData = &m_inputBuffer[frame * blockSize];
        
        // Perform FFT
        m_fft->forward(frameData, m_fftOut.data());
        
        // Compute magnitudes and accumulate
        for (size_t i = 1; i <= blockSize / 2; ++i) {
            double real = m_fftOut[2 * i];
            double imag = m_fftOut[2 * i + 1];
            double magnitude = std::sqrt(real * real + imag * imag);
            
            if (magnitude > m_globalMax) {
                m_globalMax = magnitude;
            }
            
            m_accumulatedSpectrum[i-1] += magnitude;
        }
        
        m_frameCount++;
    }
}

BIPlugin::FeatureSet BIPlugin::getRemainingFeatures()
{
    FeatureSet featureSet;

    // Process any remaining frames in the buffer
    if (!m_inputBuffer.empty()) {
        size_t remainingFrames = m_inputBuffer.size() / m_blockSize;
        if (remainingFrames > 0) {
            processBatch(remainingFrames);
        }
        m_inputBuffer.clear();
    }

    if (m_frameCount == 0) {
        // std::cerr << "BIPlugin: frameCount is 0" << std::endl;
        return featureSet;
    }

    // 1. Calculate Mean Spectrum (Magnitude)
    size_t numBins = m_numBins;
    std::vector<double> meanSpectrum(numBins);
    
    for (size_t b = 0; b < numBins; ++b) {
        meanSpectrum[b] = m_accumulatedSpectrum[b] / m_frameCount;
    }

    // 2. Convert to dB Normalized (Max 0)
    // Using same logic as soundecology for consistency unless Maad is proven different
    float maxAmp = m_globalMax;
    if (maxAmp <= 1e-9f) maxAmp = 1.0f; 

    std::vector<double> meanSpectrumDB(numBins);
    for (size_t b = 0; b < numBins; ++b) {
        double val = meanSpectrum[b] / maxAmp;
        if (val < 1e-12) val = 1e-12; 
        meanSpectrumDB[b] = 20.0 * std::log10(val);
    }

    // 3. Extract Segment and Normalize
    float df = (m_inputSampleRate / 2.0f) / numBins;
    float rows_width = 1.0f / df;
    
    size_t minBin = static_cast<size_t>(std::floor(m_minFreq * rows_width));
    size_t maxBin = static_cast<size_t>(std::ceil(m_maxFreq * rows_width));
    
    if (minBin >= numBins) minBin = numBins - 1;
    if (maxBin > numBins) maxBin = numBins;
    
    if (minBin >= maxBin) {
        Feature f;
        f.hasTimestamp = true;
        f.timestamp = Vamp::RealTime::frame2RealTime(0, m_inputSampleRate);
        f.values.push_back(0.0f);
        featureSet[0].push_back(f);
        return featureSet;
    }

    // Find min in segment
    double minVal = 1e9;
    for (size_t b = minBin; b < maxBin; ++b) {
        if (meanSpectrumDB[b] < minVal) minVal = meanSpectrumDB[b];
    }
    
    // 4. Calculate Area
    double bi = 0.0;
    for (size_t b = minBin; b < maxBin; ++b) {
        bi += (meanSpectrumDB[b] - minVal) * rows_width;
    }

    Feature feature;
    feature.hasTimestamp = true;
    feature.timestamp = Vamp::RealTime::frame2RealTime(0, m_inputSampleRate);
    feature.values.push_back(static_cast<float>(bi));
    
    featureSet[0].push_back(feature);

    return featureSet;
}
