
#include "THPlugin.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>

using std::string;
using std::vector;
using std::cerr;
using std::endl;

THPlugin::THPlugin(float inputSampleRate) :
    Plugin(inputSampleRate),
    m_blockSize(0),
    m_stepSize(0),
    m_fft(0),
    m_fftReal(0),
    m_complexInput(0),
    m_complexOutput(0)
{
}

THPlugin::~THPlugin()
{
    delete m_fft;
    delete m_fftReal;
    delete[] m_complexInput;
    delete[] m_complexOutput;
}

string THPlugin::getIdentifier() const { return "th"; }
string THPlugin::getName() const { return "Temporal Entropy (Ht)"; }
string THPlugin::getDescription() const { return "Calculates the Temporal Entropy (Ht) component of Acoustic Richness."; }
string THPlugin::getMaker() const { return "ReVAMP"; }
int THPlugin::getPluginVersion() const { return 1; }
string THPlugin::getCopyright() const { return "MIT License"; }

THPlugin::InputDomain THPlugin::getInputDomain() const { return TimeDomain; }
size_t THPlugin::getPreferredBlockSize() const { return 4096; }
size_t THPlugin::getPreferredStepSize() const { return 1024; }
size_t THPlugin::getMinChannelCount() const { return 1; }
size_t THPlugin::getMaxChannelCount() const { return 1; }

THPlugin::ParameterList THPlugin::getParameterDescriptors() const {
    return ParameterList();
}

float THPlugin::getParameter(string identifier) const {
    return 0;
}

void THPlugin::setParameter(string identifier, float value) {
}

THPlugin::ProgramList THPlugin::getPrograms() const { return ProgramList(); }
string THPlugin::getCurrentProgram() const { return ""; }
void THPlugin::selectProgram(string name) { }

THPlugin::OutputList THPlugin::getOutputDescriptors() const
{
    OutputList list;

    OutputDescriptor d;
    d.identifier = "th";
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

bool THPlugin::initialise(size_t channels, size_t stepSize, size_t blockSize)
{
    if (channels < getMinChannelCount() || channels > getMaxChannelCount()) return false;

    m_blockSize = blockSize;
    m_stepSize = stepSize;

    m_fft = new Vamp::FFTComplex(blockSize);
    m_fftReal = new Vamp::FFTReal(blockSize);
    m_complexInput = new double[blockSize * 2];
    m_complexOutput = new double[blockSize * 2];
    
    m_sumAmplitude = 0.0;
    m_sumAmplitudeLogAmplitude = 0.0;
    m_count = 0;

    return true;
}

void THPlugin::reset()
{
    m_sumAmplitude = 0.0;
    m_sumAmplitudeLogAmplitude = 0.0;
    m_count = 0;
}

THPlugin::FeatureSet THPlugin::process(const float *const *inputBuffers, Vamp::RealTime timestamp)
{
    if (m_stepSize == 0 || m_blockSize == 0) return FeatureSet();

    // 1. Prepare Input (Real)
    // We use m_complexInput as a temporary buffer for real input.
    // It has size 2*blockSize, so it fits blockSize doubles.
    for (size_t i = 0; i < m_blockSize; ++i) {
        m_complexInput[i] = (double)inputBuffers[0][i];
    }

    // 2. Forward FFT (Real -> Complex)
    // Output goes to m_complexOutput.
    // FFTReal produces N/2 + 1 complex values (interleaved).
    m_fftReal->forward(m_complexInput, m_complexOutput);

    // 3. Analytic Signal (Hilbert)
    // m_complexOutput contains:
    // [Re(0), Im(0), Re(1), Im(1), ..., Re(N/2), Im(N/2)]
    // Size is (N/2 + 1) * 2 doubles.
    
    // DC (0) - keep
    // Pos (1 to N/2 - 1) - double
    // Nyquist (N/2) - keep
    // Neg (N/2 + 1 to N - 1) - zero (implicit, we need to fill zeros)
    
    // k=1 to N/2 - 1: Double
    for (size_t k = 1; k < m_blockSize / 2; ++k) {
        m_complexOutput[k * 2] *= 2.0;
        m_complexOutput[k * 2 + 1] *= 2.0;
    }
    
    // Zero out the negative frequencies (from N/2 + 1 to N - 1)
    // The buffer m_complexOutput needs to be full size (2*N) for the inverse FFT.
    // FFTReal only wrote up to index 2*(N/2) + 1 = N+1.
    // So we zero from N+2 to 2*N - 1.
    
    size_t startZero = (m_blockSize / 2 + 1) * 2;
    size_t endZero = m_blockSize * 2;
    std::fill(m_complexOutput + startZero, m_complexOutput + endZero, 0.0);

    // 4. Inverse FFT (Complex -> Complex)
    // We use FFTComplex because the analytic signal is complex.
    m_fft->inverse(m_complexOutput, m_complexInput); // Output goes to m_complexInput

    // 5. Magnitude (Envelope) and Accumulate
    // Extract center stepSize samples to avoid edge effects
    size_t startOffset = 0;
    if (m_blockSize > m_stepSize) {
        startOffset = (m_blockSize - m_stepSize) / 2;
    }
    
    for (size_t i = 0; i < m_stepSize; ++i) {
        size_t idx = startOffset + i;
        if (idx >= m_blockSize) break;
        
        double re = m_complexInput[idx * 2];
        double im = m_complexInput[idx * 2 + 1];
        double mag = sqrt(re * re + im * im);
        
        if (mag > 1e-10) { // Avoid log(0)
            m_sumAmplitude += mag;
            m_sumAmplitudeLogAmplitude += mag * log(mag);
        }
        m_count++;
    }

    return FeatureSet();
}

double THPlugin::computeTH() const
{
    if (m_count <= 1 || m_sumAmplitude <= 0.0) return 0.0;

    double entropySum = (1.0 / m_sumAmplitude) * m_sumAmplitudeLogAmplitude - log(m_sumAmplitude);
    return -entropySum / log((double)m_count);
}

THPlugin::FeatureSet THPlugin::getRemainingFeatures()
{
    FeatureSet fs;
    
    if (m_count == 0) return fs;

    double Ht = computeTH();

    // Output Ht
    Feature fHt;
    fHt.hasTimestamp = true;
    fHt.timestamp = Vamp::RealTime::zeroTime;
    fHt.values.push_back((float)Ht);
    fs[0].push_back(fHt); // Output 0 is Ht

    return fs;
}
