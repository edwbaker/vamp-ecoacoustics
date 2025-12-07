#ifndef _ADI_PLUGIN_H_
#define _ADI_PLUGIN_H_

#include "ACIBasePlugin.h"

class ADIPlugin : public ACIBasePlugin
{
public:
    ADIPlugin(float inputSampleRate);
    virtual ~ADIPlugin();

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

    size_t getPreferredBlockSize() const;
    size_t getPreferredStepSize() const;

    FeatureSet getRemainingFeatures();

protected:
    // Parameters
    float m_minFreq;
    float m_maxFreq;
    float m_binStep; // Hz
    float m_dbThreshold;
    int m_indexType; // 0: Shannon, 1: Simpson, 2: Inverse Simpson

    // We need to store all data to normalize by global max
    // ACIBasePlugin already has m_spectralData
};

#endif
