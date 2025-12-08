#ifndef _SH_PLUGIN_H_
#define _SH_PLUGIN_H_

#include "EcoacousticSpectralPlugin.h"

class SHPlugin : public EcoacousticSpectralPlugin
{
public:
    SHPlugin(float inputSampleRate);
    virtual ~SHPlugin();

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

    // Public method to compute SH, useful for HPlugin
    double computeSH() const;

protected:
    void processBatch(size_t numFrames);

    // Accumulator for Mean Spectrum
    std::vector<double> m_accumulatedSpectrum;
};

#endif
