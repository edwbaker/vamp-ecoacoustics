
#include "AmplitudeIndexPlugin.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>

using std::string;
using std::vector;
using std::cerr;
using std::endl;

AmplitudeIndexPlugin::AmplitudeIndexPlugin(float inputSampleRate) :
    Plugin(inputSampleRate),
    m_blockSize(0),
    m_stepSize(0),
    m_bitDepth(16.0f),
    m_fft(0),
    m_complexInput(0),
    m_complexOutput(0)
{
}

AmplitudeIndexPlugin::~AmplitudeIndexPlugin()
{
    delete m_fft;
    delete[] m_complexInput;
    delete[] m_complexOutput;
}

string AmplitudeIndexPlugin::getIdentifier() const { return "amplitude-index"; }
string AmplitudeIndexPlugin::getName() const { return "Amplitude Index (M)"; }
string AmplitudeIndexPlugin::getDescription() const { return "Calculates the Amplitude Index (M) component of Acoustic Richness."; }
string AmplitudeIndexPlugin::getMaker() const { return "ReVAMP"; }
int AmplitudeIndexPlugin::getPluginVersion() const { return 1; }
string AmplitudeIndexPlugin::getCopyright() const { return "MIT License"; }

AmplitudeIndexPlugin::InputDomain AmplitudeIndexPlugin::getInputDomain() const { return TimeDomain; }
size_t AmplitudeIndexPlugin::getPreferredBlockSize() const { return 4096; }
size_t AmplitudeIndexPlugin::getPreferredStepSize() const { return 1024; }
size_t AmplitudeIndexPlugin::getMinChannelCount() const { return 1; }
size_t AmplitudeIndexPlugin::getMaxChannelCount() const { return 1; }

AmplitudeIndexPlugin::ParameterList AmplitudeIndexPlugin::getParameterDescriptors() const {
    ParameterList list;
    ParameterDescriptor d;
    d.identifier = "bitdepth";
    d.name = "Bit Depth";
    d.description = "Bit depth of the original audio file (used for M normalization)";
    d.unit = "bits";
    d.minValue = 1;
    d.maxValue = 64;
    d.defaultValue = 16;
    d.isQuantized = true;
    d.quantizeStep = 1;
    list.push_back(d);
    return list;
}

float AmplitudeIndexPlugin::getParameter(string identifier) const {
    if (identifier == "bitdepth") return m_bitDepth;
    return 0;
}

void AmplitudeIndexPlugin::setParameter(string identifier, float value) {
    if (identifier == "bitdepth") {
        m_bitDepth = value;
    }
}

AmplitudeIndexPlugin::ProgramList AmplitudeIndexPlugin::getPrograms() const { return ProgramList(); }
string AmplitudeIndexPlugin::getCurrentProgram() const { return ""; }
void AmplitudeIndexPlugin::selectProgram(string name) { }

AmplitudeIndexPlugin::OutputList AmplitudeIndexPlugin::getOutputDescriptors() const
{
    OutputList list;

    OutputDescriptor d;
    d.identifier = "m";
    d.name = "Amplitude Index (M)";
    d.description = "Median of the amplitude envelope";
    d.unit = "";
    d.hasFixedBinCount = true;
    d.binCount = 1;
    d.hasKnownExtents = false;
    d.isQuantized = false;
    d.sampleType = OutputDescriptor::VariableSampleRate; // Summary feature
    d.hasDuration = false;
    list.push_back(d);

    return list;
}

bool AmplitudeIndexPlugin::initialise(size_t channels, size_t stepSize, size_t blockSize)
{
    if (channels < getMinChannelCount() || channels > getMaxChannelCount()) return false;

    m_blockSize = blockSize;
    m_stepSize = stepSize;

    m_fft = new Vamp::FFTComplex(blockSize);
    m_complexInput = new double[blockSize * 2];
    m_complexOutput = new double[blockSize * 2];
    
    m_envelope.clear();
    m_envelope.reserve(100000); // Pre-allocate some space

    return true;
}

void AmplitudeIndexPlugin::reset()
{
    m_envelope.clear();
}

AmplitudeIndexPlugin::FeatureSet AmplitudeIndexPlugin::process(const float *const *inputBuffers, Vamp::RealTime timestamp)
{
    if (m_stepSize == 0 || m_blockSize == 0) return FeatureSet();

    // 1. Prepare Input (Real -> Complex)
    for (size_t i = 0; i < m_blockSize; ++i) {
        m_complexInput[i * 2] = (double)inputBuffers[0][i];
        m_complexInput[i * 2 + 1] = 0.0;
    }

    // 2. Forward FFT
    m_fft->forward(m_complexInput, m_complexOutput);

    // 3. Analytic Signal (Hilbert)
    // DC (0) - keep
    // Pos (1 to N/2 - 1) - double
    // Nyquist (N/2) - keep
    // Neg (N/2 + 1 to N - 1) - zero
    
    // k=1 to N/2 - 1: Double
    for (size_t k = 1; k < m_blockSize / 2; ++k) {
        m_complexOutput[k * 2] *= 2.0;
        m_complexOutput[k * 2 + 1] *= 2.0;
    }
    
    // k=N/2 + 1 to N - 1: Zero
    for (size_t k = m_blockSize / 2 + 1; k < m_blockSize; ++k) {
        m_complexOutput[k * 2] = 0.0;
        m_complexOutput[k * 2 + 1] = 0.0;
    }

    // 4. Inverse FFT
    m_fft->inverse(m_complexOutput, m_complexInput); // Output goes to m_complexInput

    // 5. Magnitude (Envelope) and Store
    // Extract center stepSize samples to avoid edge effects
    // If blockSize < stepSize, this logic fails, but we requested blockSize > stepSize
    size_t startOffset = 0;
    if (m_blockSize > m_stepSize) {
        startOffset = (m_blockSize - m_stepSize) / 2;
    }
    
    for (size_t i = 0; i < m_stepSize; ++i) {
        size_t idx = startOffset + i;
        if (idx >= m_blockSize) break;
        
        double re = m_complexInput[idx * 2];
        double im = m_complexInput[idx * 2 + 1];
        float mag = (float)sqrt(re * re + im * im);
        
        m_envelope.push_back(mag);
    }

    return FeatureSet();
}

AmplitudeIndexPlugin::FeatureSet AmplitudeIndexPlugin::getRemainingFeatures()
{
    FeatureSet fs;
    
    if (m_envelope.empty()) return fs;

    // Calculate M (Median)
    vector<float> sortedEnv = m_envelope;
    std::sort(sortedEnv.begin(), sortedEnv.end());
    
    float M = 0.0;
    size_t n = sortedEnv.size();
    if (n > 0) {
        if (n % 2 == 0) {
            M = (sortedEnv[n / 2 - 1] + sortedEnv[n / 2]) / 2.0f;
        } else {
            M = sortedEnv[n / 2];
        }
    }
    
    // Normalize M based on bit depth
    float scale = (float)pow(2.0, 1.0 - m_bitDepth);
    M *= scale;

    // Output M
    Feature fM;
    fM.hasTimestamp = true;
    fM.timestamp = Vamp::RealTime::zeroTime;
    fM.values.push_back(M);
    fs[0].push_back(fM); // Output 0 is M

    return fs;
}
