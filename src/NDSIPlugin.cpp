#include "NDSIPlugin.h"
#include <cmath>
#include <iostream>
#include <algorithm>

using std::string;
using std::vector;
using std::cerr;
using std::endl;

NDSIPlugin::NDSIPlugin(float inputSampleRate) :
    EcoacousticSpectralPlugin(inputSampleRate),
    m_anthroMin(1000.0f),
    m_anthroMax(2000.0f),
    m_bioMin(2000.0f),
    m_bioMax(8000.0f)
{
}

NDSIPlugin::~NDSIPlugin()
{
}

string NDSIPlugin::getIdentifier() const { return "ndsi"; }
string NDSIPlugin::getName() const { return "Normalized Difference Soundscape Index (NDSI)"; }
string NDSIPlugin::getDescription() const { return "Calculates the Normalized Difference Soundscape Index (NDSI) mirroring seewave::NDSI."; }
string NDSIPlugin::getMaker() const { return "ReVAMP"; }
int NDSIPlugin::getPluginVersion() const { return 1; }
string NDSIPlugin::getCopyright() const { return "MIT License"; }

NDSIPlugin::ParameterList NDSIPlugin::getParameterDescriptors() const
{
    ParameterList list;

    ParameterDescriptor d;
    d.identifier = "anthro_min";
    d.name = "Anthropophony Min Freq";
    d.description = "Minimum frequency for anthropophony (Hz)";
    d.unit = "Hz";
    d.minValue = 0;
    d.maxValue = m_inputSampleRate / 2;
    d.defaultValue = 1000;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "anthro_max";
    d.name = "Anthropophony Max Freq";
    d.description = "Maximum frequency for anthropophony (Hz)";
    d.defaultValue = 2000;
    list.push_back(d);

    d.identifier = "bio_min";
    d.name = "Biophony Min Freq";
    d.description = "Minimum frequency for biophony (Hz)";
    d.defaultValue = 2000;
    list.push_back(d);

    d.identifier = "bio_max";
    d.name = "Biophony Max Freq";
    d.description = "Maximum frequency for biophony (Hz)";
    d.defaultValue = 8000;
    list.push_back(d);

    return list;
}

float NDSIPlugin::getParameter(string identifier) const
{
    if (identifier == "anthro_min") return m_anthroMin;
    if (identifier == "anthro_max") return m_anthroMax;
    if (identifier == "bio_min") return m_bioMin;
    if (identifier == "bio_max") return m_bioMax;
    return 0;
}

void NDSIPlugin::setParameter(string identifier, float value)
{
    if (identifier == "anthro_min") m_anthroMin = value;
    if (identifier == "anthro_max") m_anthroMax = value;
    if (identifier == "bio_min") m_bioMin = value;
    if (identifier == "bio_max") m_bioMax = value;
}

NDSIPlugin::ProgramList NDSIPlugin::getPrograms() const { return ProgramList(); }
string NDSIPlugin::getCurrentProgram() const { return ""; }
void NDSIPlugin::selectProgram(string name) { }

NDSIPlugin::OutputList NDSIPlugin::getOutputDescriptors() const
{
    OutputList list;

    OutputDescriptor d;
    d.identifier = "ndsi";
    d.name = "NDSI";
    d.description = "Normalized Difference Soundscape Index";
    d.unit = "";
    d.hasFixedBinCount = true;
    d.binCount = 1;
    d.hasKnownExtents = true;
    d.minValue = -1;
    d.maxValue = 1;
    d.isQuantized = false;
    d.sampleType = OutputDescriptor::VariableSampleRate;
    d.hasDuration = false;
    list.push_back(d);

    OutputDescriptor d2;
    d2.identifier = "biophony";
    d2.name = "Biophony";
    d2.description = "Biophony Energy";
    d2.sampleType = OutputDescriptor::VariableSampleRate;
    list.push_back(d2);

    OutputDescriptor d3;
    d3.identifier = "anthropophony";
    d3.name = "Anthropophony";
    d3.description = "Anthropophony Energy";
    d3.sampleType = OutputDescriptor::VariableSampleRate;
    list.push_back(d3);

    return list;
}

bool NDSIPlugin::initialise(size_t channels, size_t stepSize, size_t blockSize)
{
    if (!EcoacousticSpectralPlugin::initialise(channels, stepSize, blockSize)) return false;

    m_accumulatedPower.assign(blockSize / 2 + 1, 0.0);
    return true;
}

void NDSIPlugin::reset()
{
    EcoacousticSpectralPlugin::reset();
    std::fill(m_accumulatedPower.begin(), m_accumulatedPower.end(), 0.0);
}

void NDSIPlugin::processBatch(size_t numFrames)
{
    if (numFrames == 0) return;

    size_t blockSize = m_blockSize;
    
    for (size_t frame = 0; frame < numFrames; ++frame) {
        double* frameData = &m_inputBuffer[frame * blockSize];
        
        m_fft->forward(frameData, m_fftOut.data());
        
        // Accumulate Power
        // DC (i=0)
        double dc = m_fftOut[0];
        m_accumulatedPower[0] += dc * dc;

        // Positive Frequencies (i=1 to N/2)
        for (size_t i = 1; i <= blockSize / 2; ++i) {
            double real = m_fftOut[2 * i];
            double imag = m_fftOut[2 * i + 1];
            double power = real * real + imag * imag;
            m_accumulatedPower[i] += power;
        }
    }
}

NDSIPlugin::FeatureSet NDSIPlugin::getRemainingFeatures()
{
    // Process any remaining frames in buffer
    size_t currentFrames = m_inputBuffer.size() / m_blockSize;
    if (currentFrames > 0) {
        processBatch(currentFrames);
        m_inputBuffer.clear();
    }

    FeatureSet fs;
    
    double anthroSum = 0.0;
    double bioSum = 0.0;
    
    double binWidth = m_inputSampleRate / m_blockSize;

    for (size_t i = 0; i < m_accumulatedPower.size(); ++i) {
        double freq = i * binWidth;
        
        // Check if bin is within Anthro range
        // We use >= min and < max to avoid double counting if ranges touch
        
        if (freq >= m_anthroMin && freq < m_anthroMax) {
            anthroSum += m_accumulatedPower[i];
        }
        
        if (freq >= m_bioMin && freq < m_bioMax) {
            bioSum += m_accumulatedPower[i];
        }
    }

    double ndsi = 0.0;
    if (bioSum + anthroSum > 0) {
        ndsi = (bioSum - anthroSum) / (bioSum + anthroSum);
    }

    Feature f;
    f.hasTimestamp = true;
    f.timestamp = Vamp::RealTime::zeroTime;
    f.values.push_back((float)ndsi);
    fs[0].push_back(f);

    Feature fB;
    fB.hasTimestamp = true;
    fB.timestamp = Vamp::RealTime::zeroTime;
    fB.values.push_back((float)bioSum);
    fs[1].push_back(fB);

    Feature fA;
    fA.hasTimestamp = true;
    fA.timestamp = Vamp::RealTime::zeroTime;
    fA.values.push_back((float)anthroSum);
    fs[2].push_back(fA);

    return fs;
}
