/*
    Ecoacoustic Vamp Plugins
    
    High-performance implementations of acoustic indices for 
    bioacoustics and soundscape ecology analysis.
    
    (C) Ed Baker 2025. Licensed under GPL (>=2).
*/

// HPlugin.cpp - Total Entropy (H = SH * TH)
// Matches seewave::H

#include "HPlugin.h"
#include "../ext/pocketfft_hdronly.h"
#include <cmath>
#include <algorithm>
#include <complex>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using std::string;
using std::vector;

HPlugin::HPlugin(float inputSampleRate) :
    Plugin(inputSampleRate),
    m_blockSize(0),
    m_stepSize(0),
    m_wl(512),
    m_spectrumCount(0),
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

string
HPlugin::getIdentifier() const
{
    return "h";
}

string
HPlugin::getName() const
{
    return "Total Entropy (H)";
}

string
HPlugin::getDescription() const
{
    return "Calculates the Total Entropy (H = SH * TH) matching seewave::H().";
}

string
HPlugin::getMaker() const
{
    return "Ecoacoustic-Vamp-Plugins";
}

int
HPlugin::getPluginVersion() const
{
    return 1;
}

string
HPlugin::getCopyright() const
{
    return "(C) Ed Baker 2025. Licensed under GPL (>=2).";
}

Vamp::Plugin::InputDomain
HPlugin::getInputDomain() const
{
    return TimeDomain;
}

size_t
HPlugin::getPreferredBlockSize() const
{
    return 512;
}

size_t
HPlugin::getPreferredStepSize() const
{
    return 512;
}

size_t
HPlugin::getMinChannelCount() const
{
    return 1;
}

size_t
HPlugin::getMaxChannelCount() const
{
    return 2;
}

HPlugin::ParameterList HPlugin::getParameterDescriptors() const {
    return ParameterList();
}

float HPlugin::getParameter(string identifier) const {
    (void)identifier;
    return 0;
}

void HPlugin::setParameter(string identifier, float value) {
    (void)identifier;
    (void)value;
}

HPlugin::ProgramList HPlugin::getPrograms() const { return ProgramList(); }
string HPlugin::getCurrentProgram() const { return ""; }
void HPlugin::selectProgram(string name) { (void)name; }

HPlugin::OutputList HPlugin::getOutputDescriptors() const
{
    OutputList list;

    OutputDescriptor d;
    d.identifier = "h";
    d.name = "Total Entropy (H)";
    d.description = "Product of Spectral Entropy (from meanspec) and Temporal Entropy";
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
    if (channels < getMinChannelCount() || channels > getMaxChannelCount()) return false;

    m_blockSize = blockSize;
    m_stepSize = stepSize;
    
    m_wl = 512;
    size_t numBins = m_wl / 2;
    m_sumSpectrum.resize(numBins, 0.0);
    m_spectrumCount = 0;
    
    // Hanning window for STFT (float for performance)
    m_window.resize(m_wl);
    for (size_t i = 0; i < m_wl; ++i) {
        m_window[i] = static_cast<float>(0.5 * (1.0 - std::cos(2.0 * M_PI * i / (m_wl - 1))));
    }
    
    // Initialize TH plugin
    if (!m_thPlugin->initialise(channels, stepSize, blockSize)) return false;
    
    // Samples are accumulated once, in the TH plugin; computeSH() reads them from there.
    m_inputBuffer.clear();

    return true;
}

void HPlugin::reset()
{
    std::fill(m_sumSpectrum.begin(), m_sumSpectrum.end(), 0.0);
    m_spectrumCount = 0;
    m_inputBuffer.clear();
    m_thPlugin->reset();
}

HPlugin::FeatureSet HPlugin::process(const float *const *inputBuffers, Vamp::RealTime timestamp)
{
    // The TH plugin accumulates the samples (as received, floats, in arrival order) for both indices; each is
    // widened to double where it is used
    m_thPlugin->process(inputBuffers, timestamp);

    return FeatureSet();
}

double HPlugin::computeSH() const
{
    // The samples TH accumulated (floats), widened to double here as they are used: the same values, in the same
    // order, as the separate double copy this plugin used to keep
    const SampleStore& samples = m_thPlugin->samples();
    const size_t received = samples.size();
    size_t n = received;

    // Remove trailing zeros (host may pad last block)
    n = samples.sizeWithoutTrailingZeros();
    // Keep some zeros back if we removed too many (file might actually end with zeros)
    size_t removed = received - n;
    if (removed > 0 && removed < m_blockSize) {
        // Likely actual padding, keep trimmed size
    } else {
        n = received;
    }
    
    if (n < m_wl) return 0.0;
    
    size_t numBins = m_wl / 2;
    std::vector<double> avgSpectrum(numBins, 0.0);
    size_t frameCount = 0;
    
    std::vector<double> frame(m_wl);
    // The plan pocketfft::r2c builds (and discards) for every call, built once for all the frames. Its exec runs
    // on the windowed frame in place; r2c would copy the frame to a work array, run this same exec and unpack the
    // result. The result is in half-complex order: r0, r1, i1, r2, i2, ...
    pocketfft::detail::pocketfft_r<double> plan(m_wl);

    std::vector<float> scratch(m_wl);
    const float* win = m_window.data();
    const double wlInv = 2.0 / static_cast<double>(m_wl);
    
    for (size_t start = 0; start + m_wl <= n; start += m_wl) {
        // Apply Hanning window
        const float* src = samples.view(start, m_wl, scratch.data());
        for (size_t i = 0; i < m_wl; ++i) {
            frame[i] = static_cast<double>(src[i]) * win[i];
        }
        
        // FFT
        plan.exec(frame.data(), 1.0, true);

        // Accumulate magnitude spectrum: bin 0 is (r0, 0), bin i is (frame[2i-1], frame[2i])
        {
            double mag = std::abs(std::complex<double>(frame[0], 0.0)) * wlInv;
            avgSpectrum[0] += mag;
        }
        for (size_t i = 1; i < numBins; ++i) {
            double mag = std::abs(std::complex<double>(frame[2 * i - 1], frame[2 * i])) * wlInv;
            avgSpectrum[i] += mag;
        }
        frameCount++;
    }
    
    if (frameCount == 0) return 0.0;
    
    // Average spectra
    for (size_t i = 0; i < numBins; ++i) {
        avgSpectrum[i] /= frameCount;
    }
    
    // Normalize to max=1
    double maxVal = 0.0;
    for (size_t i = 0; i < numBins; ++i) {
        if (avgSpectrum[i] > maxVal) maxVal = avgSpectrum[i];
    }
    if (maxVal == 0.0) return 0.0;
    
    for (size_t i = 0; i < numBins; ++i) {
        avgSpectrum[i] /= maxVal;
    }
    
    // Replace zeros with 1e-6
    for (size_t i = 0; i < numBins; ++i) {
        if (avgSpectrum[i] == 0.0) avgSpectrum[i] = 1e-6;
    }
    
    // Compute Shannon entropy - sh() on normalized spectrum
    const double epsilon = 1e-7;
    
    // Create PMF
    double sumSpec = 0.0;
    for (size_t i = 0; i < numBins; ++i) {
        sumSpec += avgSpectrum[i];
    }
    
    if (sumSpec == 0.0) return 0.0;
    
    std::vector<double> pmf(numBins);
    for (size_t i = 0; i < numBins; ++i) {
        pmf[i] = avgSpectrum[i] / sumSpec;
        if (pmf[i] == 0.0) pmf[i] = epsilon;
    }
    
    // Shannon entropy
    double entropySum = 0.0;
    for (size_t i = 0; i < numBins; ++i) {
        entropySum += pmf[i] * std::log(pmf[i]);
    }
    
    double sh = -entropySum / std::log(static_cast<double>(numBins));
    return sh;
}

HPlugin::FeatureSet HPlugin::getRemainingFeatures()
{
    FeatureSet fs;
    
    // Compute SH using meanspec approach
    double sh = computeSH();
    
    // Get TH from THPlugin
    FeatureSet thFs = m_thPlugin->getRemainingFeatures();
    double th = 0.0;
    if (!thFs.empty() && !thFs[0].empty()) {
        th = thFs[0][0].values[0];
    }
    
    double h = sh * th;

    Feature f;
    f.hasTimestamp = true;
    f.timestamp = Vamp::RealTime::zeroTime;
    f.values.push_back(static_cast<float>(h));
    fs[0].push_back(f);

    return fs;
}

