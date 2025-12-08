#ifndef _NDSI_PLUGIN_H_
#define _NDSI_PLUGIN_H_

#include "EcoacousticSpectralPlugin.h"

class NDSIPlugin : public EcoacousticSpectralPlugin
{
public:
    NDSIPlugin(float inputSampleRate);
    virtual ~NDSIPlugin();

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

    // Accumulator for Power Spectrum
    std::vector<double> m_accumulatedPower;

    // Parameters
    float m_anthroMin;
    float m_anthroMax;
    float m_bioMin;
    float m_bioMax;
};

#endif
