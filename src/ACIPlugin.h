// ACIPlugin.h - Acoustic Complexity Index
// Uses Vamp SDK FFT implementation

#ifndef _ACI_PLUGIN_H_
#define _ACI_PLUGIN_H_

#include "vamp-sdk/Plugin.h"
#include "vamp-sdk/FFT.h"

using std::string;

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

    FeatureSet process(const float *const *inputBuffers,
                       Vamp::RealTime timestamp);

    FeatureSet getRemainingFeatures();

protected:
    // Parameters matching seewave::ACI
    float m_minFreq;     // flim[1] in kHz
    float m_maxFreq;     // flim[2] in kHz
    int m_nbWindows;     // number of temporal windows
    
    // Internal state
    size_t m_blockSize;
    size_t m_stepSize;
    size_t m_channels;
    
    // FFT
    Vamp::FFTReal *m_fft;
    std::vector<double> m_fftOutput; // Buffer for FFT output
    
    // Storage for spectral data
    std::vector<std::vector<float>> m_spectralData;  // freq x time
    size_t m_frameCount;
};

#endif
