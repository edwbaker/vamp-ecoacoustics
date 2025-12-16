#ifndef _ACIACC_PLUGIN_H_
#define _ACIACC_PLUGIN_H_

#include "vamp-sdk/Plugin.h"
#include "vamp-sdk/FFT.h"
#include <vector>
#include <string>

using std::string;
using std::vector;

class ACIaccPlugin : public Vamp::Plugin
{
public:
    ACIaccPlugin(float inputSampleRate);
    virtual ~ACIaccPlugin();

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
    void processFrame(size_t channel, const float* input);

    // Internal
    size_t m_channels;
    size_t m_blockSize;
    size_t m_stepSize;
    size_t m_numBins;
    size_t m_frameCount;

    // Parameters
    float m_minFreq;
    float m_maxFreq;

    // FFT
    Vamp::FFTReal *m_fft;
    std::vector<double> m_fftInput;  // Reusable FFT input buffer
    std::vector<double> m_fftOut;
    std::vector<float> m_window;     // Float for memory efficiency

    enum WindowType {
        Hanning,
        Hamming
    };
    WindowType m_windowType;

private:
    float m_clusterSize; // in seconds
    
    // Streaming calculation state - per channel for stereo support
    size_t m_framesPerCluster;
    std::vector<size_t> m_totalSamplesReceived;  // Track actual samples per channel
    std::vector<size_t> m_lastBlockValidSamples; // Track valid (non-padded) samples in last block per channel
    std::vector<size_t> m_frameCount_ch;  // Frame count per channel
    std::vector<size_t> m_spectralWriteIdx_ch;  // Write position per channel
    
    // Per-channel spectral data (magnitude squared for deferred sqrt)
    std::vector<std::vector<float>> m_spectralData_ch;  // [channel][frame*numBins + bin]
};

#endif
