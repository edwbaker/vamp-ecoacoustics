#ifndef _ADI_ACC_PLUGIN_H_
#define _ADI_ACC_PLUGIN_H_

#include "EcoacousticSpectralPlugin.h"

class ADIaccPlugin : public EcoacousticSpectralPlugin
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
    void processBatch(size_t numFrames);

    // Parameters
    float m_dbThreshold;
    float m_freqStep; // Hz

    // Histogram-based optimization
    // We store a histogram of dB values for each band
    // This allows calculating the threshold count relative to Global Max at the end
    // without storing the full spectral data.
    std::vector<std::vector<int>> m_bandHistograms;
    bool m_bandsInitialized;
    std::vector<int> m_bandStartBins;
    std::vector<int> m_bandEndBins;
};

#endif
