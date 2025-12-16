
#ifndef _TH_PLUGIN_H_
#define _TH_PLUGIN_H_

#include "vamp-sdk/Plugin.h"
#include <vector>

using std::string;
using std::vector;

class THPlugin : public Vamp::Plugin
{
public:
    THPlugin(float inputSampleRate);
    virtual ~THPlugin();

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

    // Public method to compute TH, useful for HPlugin
    double computeTH() const;

protected:
    size_t m_blockSize;
    size_t m_stepSize;
    size_t m_totalSamples;
    
    // Accumulate all samples for full-signal processing
    std::vector<double> m_inputBuffer;
};

#endif
