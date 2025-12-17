/*
    Ecoacoustic Vamp Plugins
    
    High-performance implementations of acoustic indices for 
    bioacoustics and soundscape ecology analysis.
    
    (C) Ed Baker 2025. Licensed under GPL (>=2).
*/

#ifndef _ACI_PLUGIN_H_
#define _ACI_PLUGIN_H_

#include "vamp-sdk/Plugin.h"
#include "vamp-sdk/FFT.h"
#include <vector>
#include <string>

using std::string;
using std::vector;

class ACIPlugin : public Vamp::Plugin
{
public:
    ACIPlugin(float inputSampleRate);
    virtual ~ACIPlugin();

    string getIdentifier() const;
    string getName() const;
    string getDescription() const;
    string getMaker() const;
    int getPluginVersion() const;
    string getCopyright() const;

    InputDomain getInputDomain() const;
    size_t getPreferredBlockSize() const;
    size_t getPreferredStepSize() const;
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
    FeatureSet getRemainingFeatures();

protected:
    void processFrame(const float* input);
    int m_nbWindows;

    // Internal
    size_t m_channels;
    size_t m_blockSize;
    size_t m_stepSize;
    size_t m_numBins;
    size_t m_frameCount;
    float m_globalMaxSq; // Track global maximum magnitude squared

    // Parameters
    float m_minFreq;
    float m_maxFreq;

    // FFT
    Vamp::FFTReal *m_fft;
    std::vector<double> m_fftInput;   // Reusable FFT input buffer
    std::vector<double> m_fftOut;
    std::vector<float> m_window;

    // Streaming ACI computation (for nbWindows=1)
    std::vector<float> m_prevFrame;      // Previous frame magnitudes squared
    std::vector<float> m_sumIntensity;   // Running sum of intensities per bin
    std::vector<float> m_sumAbsDiff;     // Running sum of |diff| per bin
    bool m_hasFirstFrame;
    
    // Multi-window mode: store all spectral data
    std::vector<float> m_spectralData;
    size_t m_spectralWriteIdx;
    
    enum WindowType {
        Hanning,
        Hamming
    };
    WindowType m_windowType;
};

#endif
