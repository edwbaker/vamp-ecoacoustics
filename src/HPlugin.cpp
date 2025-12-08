#include "HPlugin.h"
#include <iostream>

using std::string;
using std::vector;
using std::cerr;
using std::endl;

HPlugin::HPlugin(float inputSampleRate) :
    SHPlugin(inputSampleRate),
    m_thPlugin(nullptr)
{
    m_thPlugin = new THPlugin(inputSampleRate);
}

HPlugin::~HPlugin()
{
    if (m_thPlugin) {
        delete m_thPlugin;
        m_thPlugin = nullptr;
    }
}

string HPlugin::getIdentifier() const { return "h"; }
string HPlugin::getName() const { return "Total Entropy (H)"; }
string HPlugin::getDescription() const { return "Calculates the Total Entropy (H = SH * TH) of a signal, matching seewave implementation."; }
string HPlugin::getMaker() const { return "ReVAMP"; }
int HPlugin::getPluginVersion() const { return 1; }
string HPlugin::getCopyright() const { return "MIT License"; }

size_t HPlugin::getPreferredBlockSize() const { return SHPlugin::getPreferredBlockSize(); }
size_t HPlugin::getPreferredStepSize() const { return SHPlugin::getPreferredStepSize(); }

HPlugin::OutputList HPlugin::getOutputDescriptors() const
{
    OutputList list;

    OutputDescriptor d;
    d.identifier = "h";
    d.name = "Total Entropy (H)";
    d.description = "Product of Spectral Entropy and Temporal Entropy";
    d.unit = "";
    d.hasFixedBinCount = true;
    d.binCount = 1;
    d.hasKnownExtents = false;
    d.isQuantized = false;
    d.sampleType = OutputDescriptor::VariableSampleRate;
    d.hasDuration = false;
    list.push_back(d);

    return list;
}

bool HPlugin::initialise(size_t channels, size_t stepSize, size_t blockSize)
{
    if (!SHPlugin::initialise(channels, stepSize, blockSize)) return false;
    if (!m_thPlugin->initialise(channels, stepSize, blockSize)) return false;
    return true;
}

void HPlugin::reset()
{
    SHPlugin::reset();
    m_thPlugin->reset();
}

HPlugin::FeatureSet HPlugin::process(const float *const *inputBuffers, Vamp::RealTime timestamp)
{
    // Process TH
    m_thPlugin->process(inputBuffers, timestamp);
    
    // Process SH (base class)
    return SHPlugin::process(inputBuffers, timestamp);
}

HPlugin::FeatureSet HPlugin::getRemainingFeatures()
{
    FeatureSet fs;
    
    // Get SH
    // SHPlugin::getRemainingFeatures() flushes the buffer and returns SH in a FeatureSet.
    FeatureSet shFs = SHPlugin::getRemainingFeatures();
    double sh = 0.0;
    if (!shFs.empty() && !shFs[0].empty()) {
        sh = shFs[0][0].values[0];
    }
    
    // Get TH
    double th = m_thPlugin->computeTH();
    
    double h = sh * th;

    Feature f;
    f.hasTimestamp = true;
    f.timestamp = Vamp::RealTime::zeroTime;
    f.values.push_back((float)h);
    fs[0].push_back(f);

    return fs;
}
