/*
    Ecoacoustic Vamp Plugins
    
    High-performance implementations of acoustic indices for 
    bioacoustics and soundscape ecology analysis.
    
    (C) Ed Baker 2025. Licensed under GPL (>=2).
*/

#ifndef _H_PLUGIN_H_
#define _H_PLUGIN_H_

#include "vamp-sdk/Plugin.h"
#include "THPlugin.h"
#include <vector>

using std::string;

// HPlugin computes Total Entropy H = SH * TH matching seewave::H()
// This uses meanspec() for SH (STFT averaged spectrum) not spec() (full FFT)
class HPlugin : public Vamp::Plugin
{
public:
    HPlugin(float inputSampleRate);
    virtual ~HPlugin();

    string getIdentifier() const;
    string getName() const;
    string getDescription() const;
    string getMaker() const;
    int getPluginVersion() const;
    string getCopyright() const;

    InputDomain getInputDomain() const;
    size_t getPreferredBlockSize() const;
    size_t getPreferredStepSize() const;
    size_t getMinChannelCount() const;
    size_t getMaxChannelCount() const;

    ParameterList getParameterDescriptors() const;
    float getParameter(string identifier) const;
    void setParameter(string identifier, float value);

    ProgramList getPrograms() const;
    string getCurrentProgram() const;
    void selectProgram(string name);

    OutputList getOutputDescriptors() const;

    bool initialise(size_t channels, size_t stepSize, size_t blockSize);
    void reset();

    FeatureSet process(const float *const *inputBuffers,
                       Vamp::RealTime timestamp);

    FeatureSet getRemainingFeatures();

protected:
    size_t m_blockSize;
    size_t m_stepSize;
    size_t m_wl;  // Window length for meanspec (default 512)
    
    // Accumulated spectra for meanspec
    std::vector<double> m_sumSpectrum;
    size_t m_spectrumCount;
    std::vector<float> m_window;
    
    // For TH computation
    THPlugin* m_thPlugin;
    
    // Buffer for accumulating samples for TH
    std::vector<double> m_inputBuffer;
    
    double computeSH() const;
};

#endif
