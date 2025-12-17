/*
    Ecoacoustic Vamp Plugins
    
    High-performance implementations of acoustic indices for 
    bioacoustics and soundscape ecology analysis.
    
    (C) Ed Baker 2025. Licensed under GPL (>=2).
*/

#ifndef _ADI_ACC_PLUGIN_H_
#define _ADI_ACC_PLUGIN_H_

#include "vamp-sdk/Plugin.h"
#include "vamp-sdk/FFT.h"
#include <vector>
#include <string>

using std::string;
using std::vector;

class ADIaccPlugin : public Vamp::Plugin
{
public:
    ADIaccPlugin(float inputSampleRate);
    virtual ~ADIaccPlugin();

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
    std::vector<size_t> m_sampleCount_ch; // Track actual samples received per channel

    // FFT
    Vamp::FFTReal *m_fft;
    std::vector<double> m_fftOut;
    std::vector<double> m_fftInput; 
    std::vector<float> m_window;
    std::vector<std::vector<int>> m_binToBands;  // Pre-computed bin-to-bands mapping (each bin can be in multiple bands)

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
    float m_dbThreshold;
    float m_freqStep; // Hz

    // Histogram-based optimization - per channel
    std::vector<std::vector<std::vector<int>>> m_bandHistograms_ch;
    bool m_bandsInitialized;
    std::vector<int> m_bandStartBins;
    std::vector<int> m_bandEndBins;
};

#endif
