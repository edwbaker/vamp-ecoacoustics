#ifndef _AEI_ACC_PLUGIN_H_
#define _AEI_ACC_PLUGIN_H_

#include "ACIBasePlugin.h"

class AEIaccPlugin : public ACIBasePlugin
{
public:
    AEIaccPlugin(float inputSampleRate);
    virtual ~AEIaccPlugin();

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
};

#endif
