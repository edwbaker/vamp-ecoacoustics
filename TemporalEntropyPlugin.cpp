
#include "TemporalEntropyPlugin.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>

using std::string;
using std::vector;
using std::cerr;
using std::endl;

TemporalEntropyPlugin::TemporalEntropyPlugin(float inputSampleRate) :
    Plugin(inputSampleRate),
    m_blockSize(0),
    m_stepSize(0),
    m_fft(0),
    m_complexInput(0),
    m_complexOutput(0)
{
}

TemporalEntropyPlugin::~TemporalEntropyPlugin()
{
    delete m_fft;
    delete[] m_complexInput;
    delete[] m_complexOutput;
}

string TemporalEntropyPlugin::getIdentifier() const { return "temporal-entropy"; }
string TemporalEntropyPlugin::getName() const { return "Temporal Entropy (Ht)"; }
string TemporalEntropyPlugin::getDescription() const { return "Calculates the Temporal Entropy (Ht) component of Acoustic Richness."; }
string TemporalEntropyPlugin::getMaker() const { return "ReVAMP"; }
int TemporalEntropyPlugin::getPluginVersion() const { return 1; }
string TemporalEntropyPlugin::getCopyright() const { return "MIT License"; }

TemporalEntropyPlugin::InputDomain TemporalEntropyPlugin::getInputDomain() const { return TimeDomain; }
size_t TemporalEntropyPlugin::getPreferredBlockSize() const { return 4096; }
size_t TemporalEntropyPlugin::getPreferredStepSize() const { return 1024; }
size_t TemporalEntropyPlugin::getMinChannelCount() const { return 1; }
size_t TemporalEntropyPlugin::getMaxChannelCount() const { return 1; }

TemporalEntropyPlugin::ParameterList TemporalEntropyPlugin::getParameterDescriptors() const {
    return ParameterList();
}

float TemporalEntropyPlugin::getParameter(string identifier) const {
    return 0;
}

void TemporalEntropyPlugin::setParameter(string identifier, float value) {
}

TemporalEntropyPlugin::ProgramList TemporalEntropyPlugin::getPrograms() const { return ProgramList(); }
string TemporalEntropyPlugin::getCurrentProgram() const { return ""; }
void TemporalEntropyPlugin::selectProgram(string name) { }

TemporalEntropyPlugin::OutputList TemporalEntropyPlugin::getOutputDescriptors() const
{
    OutputList list;

    OutputDescriptor d;
    d.identifier = "ht";
    d.name = "Temporal Entropy (Ht)";
    d.description = "Temporal entropy of the amplitude envelope";
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

bool TemporalEntropyPlugin::initialise(size_t channels, size_t stepSize, size_t blockSize)
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

void TemporalEntropyPlugin::reset()
{
    m_envelope.clear();
}

TemporalEntropyPlugin::FeatureSet TemporalEntropyPlugin::process(const float *const *inputBuffers, Vamp::RealTime timestamp)
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

TemporalEntropyPlugin::FeatureSet TemporalEntropyPlugin::getRemainingFeatures()
{
    FeatureSet fs;
    
    if (m_envelope.empty()) return fs;

    // Calculate Ht (Temporal Entropy)
    // Ht = - sum(p * log2(p)) / log2(N)
    // p = A / sum(A)
    
    double sumA = 0.0;
    for (float v : m_envelope) {
        sumA += v;
    }
    
    double Ht = 0.0;
    size_t n = m_envelope.size();
    
    if (sumA > 0.0) {
        double entropySum = 0.0;
        for (float v : m_envelope) {
            if (v > 0.0f) {
                double p = v / sumA;
                entropySum += p * (log(p) / log(2.0));
            }
        }
        
        if (n > 1) {
            Ht = -entropySum / (log((double)n) / log(2.0));
        }
    }

    // Output Ht
    Feature fHt;
    fHt.hasTimestamp = true;
    fHt.timestamp = Vamp::RealTime::zeroTime;
    fHt.values.push_back((float)Ht);
    fs[0].push_back(fHt); // Output 0 is Ht

    return fs;
}
