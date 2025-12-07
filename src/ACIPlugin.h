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
    // Parameters
    float m_minFreq;
    float m_maxFreq;
    int m_nbWindows;

    // Internal
    size_t m_channels;
    size_t m_blockSize;
    size_t m_stepSize;
    int m_frameCount;

    // FFT
    Vamp::FFTReal *m_fft;
    vector<double> m_fftOut; // Buffer for FFT output (interleaved complex)

    // Batch processing
    size_t m_batchSize;
    vector<double> m_inputBuffer; // Flattened buffer for multiple frames
    vector<double> m_window; // Pre-computed window
    void processBatch(size_t numFrames);

    // Data storage for ACI calculation
    // Flattened vector storing magnitudes: Frame 0 [Bin0...BinN], Frame 1 [Bin0...BinN]
    vector<float> m_spectralData;
    size_t m_numBins; // Number of bins per frame
};

#endif
