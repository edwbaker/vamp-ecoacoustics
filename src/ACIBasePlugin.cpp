#include "ACIBasePlugin.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

ACIBasePlugin::ACIBasePlugin(float inputSampleRate) :
    Plugin(inputSampleRate),
    m_minFreq(0),
    m_maxFreq(0),
    m_channels(0),
    m_blockSize(0),
    m_stepSize(0),
    m_numBins(0),
    m_frameCount(0),
    m_fft(nullptr),
    m_batchSize(256)
{
}

ACIBasePlugin::~ACIBasePlugin()
{
    if (m_fft) {
        delete m_fft;
        m_fft = nullptr;
    }
}

Vamp::Plugin::InputDomain
ACIBasePlugin::getInputDomain() const
{
    return TimeDomain;
}

size_t
ACIBasePlugin::getPreferredBlockSize() const
{
    return 512;
}

size_t
ACIBasePlugin::getPreferredStepSize() const
{
    return 512;
}

size_t
ACIBasePlugin::getMinChannelCount() const
{
    return 1;
}

size_t
ACIBasePlugin::getMaxChannelCount() const
{
    return 1;
}

bool
ACIBasePlugin::initialise(size_t channels, size_t stepSize, size_t blockSize)
{
    if (channels < getMinChannelCount() ||
        channels > getMaxChannelCount()) return false;

    m_channels = channels;
    m_blockSize = blockSize;
    m_stepSize = stepSize;
    m_frameCount = 0;
    m_spectralData.clear();
    m_numBins = blockSize / 2; // We store bins 1 to N/2 (skipping DC)

    // Initialize FFT
    if (m_fft) delete m_fft;
    m_fft = new Vamp::FFTReal(blockSize);
    m_fftOut.resize(blockSize + 2);
    
    // Reserve input buffer
    m_inputBuffer.reserve(blockSize * m_batchSize);
    
    // Pre-compute Hamming window
    m_window.resize(m_blockSize);
    for (size_t i = 0; i < m_blockSize; ++i) {
        m_window[i] = 0.54 - 0.46 * std::cos(2.0 * M_PI * i / (m_blockSize - 1));
    }

    return true;
}

void
ACIBasePlugin::reset()
{
    m_spectralData.clear();
    m_inputBuffer.clear();
    m_frameCount = 0;
}

Vamp::Plugin::FeatureSet
ACIBasePlugin::process(const float *const *inputBuffers, Vamp::RealTime timestamp)
{
    FeatureSet fs;
    
    // Apply window and append to input buffer
    for (size_t i = 0; i < m_blockSize; ++i) {
        m_inputBuffer.push_back(inputBuffers[0][i] * m_window[i]);
    }
    
    // Check if we have enough frames for a batch
    size_t currentFrames = m_inputBuffer.size() / m_blockSize;
    if (currentFrames >= m_batchSize) {
        processBatch(currentFrames);
        m_inputBuffer.clear();
    }
    
    return fs;
}

void ACIBasePlugin::computeMagnitudes(size_t numFrames)
{
    size_t blockSize = m_blockSize;
    
    // Ensure capacity
    size_t requiredSize = m_spectralData.size() + numFrames * m_numBins;
    if (m_spectralData.capacity() < requiredSize) {
        m_spectralData.reserve(std::max(requiredSize, m_spectralData.capacity() * 2));
    }

    for (size_t frame = 0; frame < numFrames; ++frame) {
        // Get pointer to current frame in input buffer
        double* frameData = &m_inputBuffer[frame * blockSize];
        
        // Perform FFT
        m_fft->forward(frameData, m_fftOut.data());
        
        // Compute magnitudes and append to flattened vector
        // m_fftOut contains interleaved complex data: real, imag, real, imag...
        // Skip DC (i=0), start from i=1 up to Nyquist (i=blockSize/2)
        for (size_t i = 1; i <= blockSize / 2; ++i) {
            double real = m_fftOut[2 * i];
            double imag = m_fftOut[2 * i + 1];
            // Use float for storage to save memory/bandwidth
            float magnitude = static_cast<float>(std::sqrt(real * real + imag * imag));
            m_spectralData.push_back(magnitude);
        }
        
        m_frameCount++;
    }
}

void ACIBasePlugin::processBatch(size_t numFrames)
{
    if (numFrames == 0) return;
    computeMagnitudes(numFrames);
}
