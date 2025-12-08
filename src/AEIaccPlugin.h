#ifndef _AEI_ACC_PLUGIN_H_
#define _AEI_ACC_PLUGIN_H_

#include "EcoacousticSpectralPlugin.h"

class AEIaccPlugin : public EcoacousticSpectralPlugin
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
    void processBatch(size_t numFrames);

    // Parameters
    float m_dbThreshold;
    float m_freqStep; // Hz

    // Accumulators
    // Histogram: [band][db_bin]
    // dB range: -120 to 0. Resolution: 0.1 dB -> 1200 bins.
    // Bin 0: -120dB, Bin 1200: 0dB.
    // Values < -120dB go to bin 0.
    std::vector<std::vector<int>> m_bandHistograms;
    std::vector<int> m_bandStartBins;
    std::vector<int> m_bandEndBins;
    bool m_bandsInitialized;
};

#endif
