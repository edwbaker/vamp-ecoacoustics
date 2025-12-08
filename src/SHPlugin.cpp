// SHPlugin.cpp - Spectral Entropy
// Matches seewave::SH behavior

#include "SHPlugin.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

SHPlugin::SHPlugin(float inputSampleRate) :
    EcoacousticSpectralPlugin(inputSampleRate)
{
}

SHPlugin::~SHPlugin()
{
}

string
SHPlugin::getIdentifier() const
{
    return "sh";
}

string
SHPlugin::getName() const
{
    return "Spectral Entropy (seewave)";
}

string
SHPlugin::getDescription() const
{
    return "Calculates the Spectral Entropy (SH) of a signal, matching seewave implementation.";
}

string
SHPlugin::getMaker() const
{
    return "ReVAMP";
}

int
SHPlugin::getPluginVersion() const
{
    return 1;
}

string
SHPlugin::getCopyright() const
{
    return "MIT License";
}

Vamp::Plugin::ParameterList
SHPlugin::getParameterDescriptors() const
{
    ParameterList list;
    return list;
}

float
SHPlugin::getParameter(string identifier) const
{
    return 0;
}

void
SHPlugin::setParameter(string identifier, float value)
{
}

Vamp::Plugin::ProgramList
SHPlugin::getPrograms() const
{
    ProgramList list;
    return list;
}

string
SHPlugin::getCurrentProgram() const
{
    return "";
}

void
SHPlugin::selectProgram(string name)
{
}

Vamp::Plugin::OutputList
SHPlugin::getOutputDescriptors() const
{
    OutputList list;

    OutputDescriptor d;
    d.identifier = "sh";
    d.name = "Spectral Entropy";
    d.description = "Spectral Entropy (Shannon entropy of the spectral distribution)";
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
SHPlugin::getPreferredBlockSize() const
{
    return 512;
}

size_t
SHPlugin::getPreferredStepSize() const
{
    return 512; // No overlap
}

bool
SHPlugin::initialise(size_t channels, size_t stepSize, size_t blockSize)
{
    if (!EcoacousticSpectralPlugin::initialise(channels, stepSize, blockSize)) return false;

    // We include DC component to match seewave::sh behavior
    m_numBins = blockSize / 2 + 1;
    m_accumulatedSpectrum.assign(m_numBins, 0.0);

    return true;
}

void SHPlugin::processBatch(size_t numFrames)
{
    if (numFrames == 0) return;

    size_t blockSize = m_blockSize;
    
    for (size_t frame = 0; frame < numFrames; ++frame) {
        double* frameData = &m_inputBuffer[frame * blockSize];
        m_fft->forward(frameData, m_fftOut.data());
        
        // Compute magnitudes and accumulate
        // Include DC (i=0)
        for (size_t i = 0; i <= blockSize / 2; ++i) {
            double real = m_fftOut[2 * i];
            double imag = m_fftOut[2 * i + 1];
            double magnitude = std::sqrt(real * real + imag * imag);
            
            m_accumulatedSpectrum[i] += magnitude;
        }
        
        m_frameCount++;
    }
}

Vamp::Plugin::FeatureSet
SHPlugin::getRemainingFeatures()
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

    double sh = computeSH();
    
    Feature f;
    f.hasTimestamp = true;
    f.timestamp = Vamp::RealTime::frame2RealTime(0, m_inputSampleRate);
    f.values.push_back(static_cast<float>(sh));
    
    fs[0].push_back(f);
    
    return fs;
}

double SHPlugin::computeSH() const
{
    if (m_frameCount == 0) return 0.0;

    // 1. Calculate Mean Spectrum
    // seewave::sh uses the mean spectrum of the STFT
    size_t numBins = m_numBins;
    // We don't want to allocate a huge vector here if we can avoid it, 
    // but we need to normalize.
    // Actually we can compute sumSpectrum first.
    
    double sumSpectrum = 0.0;
    for (size_t b = 0; b < numBins; ++b) {
        sumSpectrum += m_accumulatedSpectrum[b]; // Sum of accumulated is proportional to sum of means
    }
    // sum(mean) = sum(acc / N) = sum(acc) / N
    
    double sh = 0.0;
    
    if (sumSpectrum > 0.0) {
        double entropySum = 0.0;
        double invFrameCount = 1.0 / m_frameCount;
        double invSumSpectrum = 1.0 / (sumSpectrum * invFrameCount); // 1 / sum(mean)

        for (size_t b = 0; b < numBins; ++b) {
            if (m_accumulatedSpectrum[b] > 0.0) {
                double meanVal = m_accumulatedSpectrum[b] * invFrameCount;
                double p = meanVal * invSumSpectrum;
                entropySum += p * (std::log(p) / std::log(2.0));
            }
        }
        
        if (numBins > 1) {
            sh = -entropySum / (std::log((double)numBins) / std::log(2.0));
        }
    }
    return sh;
}

