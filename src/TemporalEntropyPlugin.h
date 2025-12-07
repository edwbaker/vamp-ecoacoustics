
#ifndef _TEMPORAL_ENTROPY_PLUGIN_H_
#define _TEMPORAL_ENTROPY_PLUGIN_H_

#include "vamp-sdk/Plugin.h"
#include "vamp-sdk/FFT.h"
#include <vector>

using std::string;
using std::vector;

class TemporalEntropyPlugin : public Vamp::Plugin
{
public:
    TemporalEntropyPlugin(float inputSampleRate);
    virtual ~TemporalEntropyPlugin();

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
    size_t m_blockSize;
    size_t m_stepSize;
    
    // Buffer for accumulating envelope values
    vector<float> m_envelope;
    
    // FFT object
    Vamp::FFTComplex *m_fft;
    
    // Buffers for FFT
    double *m_complexInput;
    double *m_complexOutput;
};

#endif
