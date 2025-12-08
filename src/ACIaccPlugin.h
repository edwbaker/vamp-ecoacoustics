#ifndef _ACIACC_PLUGIN_H_
#define _ACIACC_PLUGIN_H_

#include "EcoacousticSpectralPlugin.h"

class ACIaccPlugin : public EcoacousticSpectralPlugin
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

    ParameterList getParameterDescriptors() const;
    float getParameter(string identifier) const;
    void setParameter(string identifier, float value);

    ProgramList getPrograms() const;
    string getCurrentProgram() const;
    void selectProgram(string name);

    OutputList getOutputDescriptors() const;

    bool initialise(size_t channels, size_t stepSize, size_t blockSize);
    void reset();

    FeatureSet getRemainingFeatures();

protected:
    void processBatch(size_t numFrames);

private:
    float m_clusterSize; // in seconds
    
    // Streaming calculation state
    float m_runningTotalACI;
    size_t m_framesPerCluster;
    
    // Reusable accumulators to avoid allocation in loop
    std::vector<float> m_sumIntensity;
    std::vector<float> m_sumAbsDiff;
};

#endif
