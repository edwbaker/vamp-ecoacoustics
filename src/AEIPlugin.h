#ifndef _AEI_PLUGIN_H_
#define _AEI_PLUGIN_H_

#include "EcoacousticSpectralPlugin.h"

class AEIPlugin : public EcoacousticSpectralPlugin
{
public:
    AEIPlugin(float inputSampleRate);
    virtual ~AEIPlugin();

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
    float m_binStep; // Hz
    float m_dbThreshold;

    // Histogram-based optimization
    std::vector<std::vector<int>> m_bandHistograms;
    bool m_bandsInitialized;
    std::vector<int> m_bandStartBins;
    std::vector<int> m_bandEndBins;
};

#endif
