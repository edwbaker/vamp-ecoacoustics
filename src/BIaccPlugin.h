/*
    Ecoacoustic Vamp Plugins
    
    High-performance implementations of acoustic indices for 
    bioacoustics and soundscape ecology analysis.
    
    (C) Ed Baker 2025. Licensed under GPL (>=2).
*/

#ifndef _BI_ACC_PLUGIN_H_
#define _BI_ACC_PLUGIN_H_

#include "vamp-sdk/Plugin.h"
#include "vamp-sdk/FFT.h"
#include <vector>
#include <string>

using std::string;
using std::vector;

class BIaccPlugin : public Vamp::Plugin
{
public:
    BIaccPlugin(float inputSampleRate);
    virtual ~BIaccPlugin();

    string getIdentifier() const;
    string getName() const;
    string getDescription() const;
    string getMaker() const;
    int getPluginVersion() const;
    string getCopyright() const;

    InputDomain getInputDomain() const;
    size_t getMinChannelCount() const;
    size_t getMaxChannelCount() const;

    ParameterList getParameterDescriptors() const;
    float getParameter(string identifier) const;
    void setParameter(string identifier, float value);

    ProgramList getPrograms() const;
    string getCurrentProgram() const;
    void selectProgram(string name);

    OutputList getOutputDescriptors() const;

    bool initialise(size_t channels, size_t stepSize, size_t blockSize);
    void reset();
    FeatureSet process(const float *const *inputBuffers, Vamp::RealTime timestamp);

    size_t getPreferredBlockSize() const;
    size_t getPreferredStepSize() const;

    FeatureSet getRemainingFeatures();

protected:
    void processBatch(size_t channel, size_t numFrames);

    // Internal
    size_t m_channels;
    size_t m_blockSize;
    size_t m_stepSize;
    size_t m_numBins;
    
    // Per-channel tracking
    std::vector<size_t> m_frameCount_ch;
    std::vector<float> m_globalMax_ch; // Track global maximum magnitude per channel

    // FFT
    Vamp::FFTReal *m_fft;
    std::vector<double> m_fftOut;
    std::vector<float> m_window;

    // Batch processing - per channel
    size_t m_batchSize;
    std::vector<std::vector<double>> m_inputBuffer_ch; 
    std::vector<std::vector<float>> m_spectralData_ch;
    
    enum WindowType {
        Hanning,
        Hamming
    };
    WindowType m_windowType;

    // Parameters
    float m_minFreq;
    float m_maxFreq;
    
    // Pre-computed frequency bin range
    size_t m_minBin;
    size_t m_maxBin;
    float m_rowsWidth;
    
    // Accumulator - per channel
    std::vector<std::vector<double>> m_accumulatedSpectrum_ch;
};

#endif
