#ifndef _BI_PLUGIN_H_
#define _BI_PLUGIN_H_

#include "EcoacousticSpectralPlugin.h"

class BIPlugin : public EcoacousticSpectralPlugin
{
public:
    BIPlugin(float inputSampleRate);
    virtual ~BIPlugin();

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
    void processBatch(size_t numFrames);

    // Parameters
    float m_minFreq;
    float m_maxFreq;

    // Accumulator
    std::vector<double> m_accumulatedSpectrum;
};

#endif
