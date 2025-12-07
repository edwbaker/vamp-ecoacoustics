#ifndef _ACI_BASE_PLUGIN_H_
#define _ACI_BASE_PLUGIN_H_

#include "vamp-sdk/Plugin.h"
#include "vamp-sdk/FFT.h"
#include <vector>
#include <string>

using std::string;
using std::vector;

class ACIBasePlugin : public Vamp::Plugin
{
public:
    ACIBasePlugin(float inputSampleRate);
    virtual ~ACIBasePlugin();

    // Common Plugin Interface
    InputDomain getInputDomain() const;
    size_t getPreferredBlockSize() const;
    size_t getPreferredStepSize() const;
    size_t getMinChannelCount() const;
    size_t getMaxChannelCount() const;

    bool initialise(size_t channels, size_t stepSize, size_t blockSize);
    void reset();

    FeatureSet process(const float *const *inputBuffers, Vamp::RealTime timestamp);

protected:
    // Helper to run FFT on buffered data and append to m_spectralData
    void computeMagnitudes(size_t numFrames);
    
    // Abstract method to be implemented by subclasses to handle batch processing
    virtual void processBatch(size_t numFrames);

    // Parameters
    float m_minFreq;
    float m_maxFreq;

    // Internal
    size_t m_channels;
    size_t m_blockSize;
    size_t m_stepSize;
    size_t m_numBins;
    size_t m_frameCount;
    float m_globalMax; // Track global maximum magnitude

    // FFT
    Vamp::FFTReal *m_fft;
    vector<double> m_fftOut;
    vector<double> m_window;

    // Batch processing
    size_t m_batchSize;
    vector<double> m_inputBuffer; 
    
    // Data storage
    vector<float> m_spectralData;
};

#endif
