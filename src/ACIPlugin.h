#ifndef _ACI_PLUGIN_H_
#define _ACI_PLUGIN_H_

#include "ACIBasePlugin.h"

class ACIPlugin : public ACIBasePlugin
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

    ParameterList getParameterDescriptors() const;
    float getParameter(string identifier) const;
    void setParameter(string identifier, float value);

    ProgramList getPrograms() const;
    string getCurrentProgram() const;
    void selectProgram(string name);

    OutputList getOutputDescriptors() const;

    bool initialise(size_t channels, size_t stepSize, size_t blockSize);

    FeatureSet getRemainingFeatures();

protected:
    int m_nbWindows;
};

#endif
