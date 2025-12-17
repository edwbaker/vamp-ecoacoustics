#ifndef _NDSI_PLUGIN_H_
#define _NDSI_PLUGIN_H_

#include "vamp-sdk/Plugin.h"
#include "vamp-sdk/FFT.h"
#include <vector>
#include <string>

using std::string;
using std::vector;

class NDSIPlugin : public Vamp::Plugin
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

    InputDomain getInputDomain() const;
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
    FeatureSet process(const float *const *inputBuffers, Vamp::RealTime timestamp);

    size_t getPreferredBlockSize() const;
    size_t getPreferredStepSize() const;

    FeatureSet getRemainingFeatures();

protected:
    void processOneSecond(const vector<double>& samples, size_t channel);
    void calculateNDSI(size_t channel, double& ndsi, double& biophony, double& anthrophony);

    size_t m_channels;
    size_t m_blockSize;
    size_t m_stepSize;
    size_t m_fftSize;      // Fixed at 1024 for pwelch
    size_t m_numBins;      // fftSize/2 = 512 (matching R's pwelch output)
    vector<size_t> m_secondCount;  // Number of 1-second chunks processed per channel
    
    float m_anthroMin;
    float m_anthroMax;
    float m_bioMin;
    float m_bioMax;

    Vamp::FFTReal *m_fft;
    vector<double> m_fftOut;
    vector<float> m_window;       // Hamming window
    double m_windowSumSq;
    
    vector<double> m_timeBuf;
    vector<double> m_detrendBuf;
    vector<double> m_windowPSD;
    
    // Pre-computed detrending constants
    double m_detrendSumT;
    double m_detrendSumTT;
    double m_detrendDenom;
    double m_detrendNfftInv;
    
    // Pre-computed PSD normalization
    double m_psdNorm;
    double m_psdNyquistNorm;
    
    vector<vector<double>> m_sampleBuffer; // Buffer to collect samples per channel
    vector<vector<double>> m_accumulatedPSD; // Sum of PSDs from each 1-second chunk per channel
    
    // Padding detection
    size_t m_totalSamplesReceived;
    size_t m_lastBlockValidSamples;
};

#endif
