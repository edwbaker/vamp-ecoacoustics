#ifndef _ADI_ACC_PLUGIN_H_
#define _ADI_ACC_PLUGIN_H_

#include "ACIBasePlugin.h"

class ADIaccPlugin : public ACIBasePlugin
{
public:
    ADIaccPlugin(float inputSampleRate);
    virtual ~ADIaccPlugin();

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
    float m_dbThreshold;
    float m_freqStep; // Hz

    // We need to store all data to normalize by global max (if required)
    // ACIBasePlugin already has m_spectralData
};

#endif
