#ifndef _SH_PLUGIN_H_
#define _SH_PLUGIN_H_

#include "vamp-sdk/Plugin.h"
#include "vamp-sdk/FFT.h"
#include <vector>
#include <string>

using std::string;
using std::vector;

class SHPlugin : public Vamp::Plugin
{
public:
    SHPlugin(float inputSampleRate);
    virtual ~SHPlugin();

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

    // Public method to compute SH, useful for HPlugin
    double computeSH() const;

protected:
    void processBatch(size_t numFrames);

    // Internal
    size_t m_channels;
    size_t m_blockSize;
    size_t m_stepSize;
    size_t m_numBins;
    size_t m_frameCount;
    float m_globalMax; // Track global maximum magnitude

    // FFT
    Vamp::FFTReal *m_fft;
    std::vector<double> m_fftOut;
    std::vector<double> m_window;

    // Batch processing
    size_t m_batchSize;
    std::vector<double> m_inputBuffer; 
    std::vector<float> m_spectralData;
    
    enum WindowType {
        Hanning,
        Hamming
    };
    WindowType m_windowType;

    // Parameters
    float m_minFreq;
    float m_maxFreq;

    // Accumulator for Mean Spectrum
    std::vector<double> m_accumulatedSpectrum;
    
    // Track actual sample count (excluding padding)
    size_t m_totalSamples;
};

#endif
