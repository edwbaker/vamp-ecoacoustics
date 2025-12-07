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

private:
    void processBatch(size_t numFrames);

    float m_minFreq;
    float m_maxFreq;
    float m_clusterSize; // in seconds

    size_t m_channels;
    size_t m_blockSize;
    size_t m_stepSize;
    size_t m_numBins;
    
    // FFT
    Vamp::FFTReal* m_fft;
    std::vector<double> m_fftOut;
    std::vector<double> m_window;
    
    // Data storage
    std::vector<double> m_inputBuffer; // Stores time-domain samples
    std::vector<float> m_spectralData; // Flattened spectral magnitudes
    
    size_t m_batchSize;
    size_t m_frameCount;
};

#endif
