/*
    Ecoacoustic Vamp Plugins
    
    High-performance implementations of acoustic indices for 
    bioacoustics and soundscape ecology analysis.
    
    (C) Ed Baker 2025. Licensed under GPL (>=2).
*/

#ifndef _AEI_ACC_PLUGIN_H_
#define _AEI_ACC_PLUGIN_H_

#include "vamp-sdk/Plugin.h"
#include "RowStore.h"
#include <vector>
#include <string>

class InPlaceRealFFT;

using std::string;
using std::vector;

class AEIaccPlugin : public Vamp::Plugin
{
public:
    AEIaccPlugin(float inputSampleRate);
    virtual ~AEIaccPlugin();

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
    void initBands();
    void processBlock(size_t channel, const float* input);

    // Internal
    size_t m_channels;
    size_t m_blockSize;
    size_t m_stepSize;
    size_t m_numBins;
    
    // Per-channel tracking
    std::vector<size_t> m_frameCount_ch;
    std::vector<double> m_globalMax_ch;      // Double: a float maximum moves the threshold by up to ~1e-6 dB,
                                             // enough to change the count for a cell that close to it
    std::vector<double> m_maxBeforeLast_ch;  // Global maximum before the last block
    std::vector<char> m_lastPadded_ch;     // The last block ended in zeros: the host padded an incomplete window

    // FFT
    InPlaceRealFFT *m_fft;
    std::vector<double> m_window;

    // Batch processing - per channel
    // Every block is transformed as it arrives; its scaled magnitudes (one row of m_numBins doubles per block)
    // are kept for the threshold counts in getRemainingFeatures().
    std::vector<RowStoreD> m_magnitude_ch;
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
    float m_freqStep;

    // Accumulators - per channel
    std::vector<std::vector<long>> m_bandAboveThresholdCount_ch;
    std::vector<std::vector<long>> m_bandTotalCount_ch;
    std::vector<int> m_bandStartBins;
    std::vector<int> m_bandEndBins;
    bool m_bandsInitialized;
};

#endif
