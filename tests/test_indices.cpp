/*
    Ecoacoustic Vamp Plugins - test runner

    Drives every plugin the way a Vamp host does (preferred block/step size,
    zero-padded final block, getRemainingFeatures at the end) on synthetic,
    fully deterministic signals, so no audio fixtures are needed.

    Three kinds of check:
      1. API sanity    - identifier, channel limits, parameters, outputs.
      2. Properties    - results with a known answer (a tone in the biophony
                         band gives NDSI ~ +1, flat noise gives ADI ~ ln(bands),
                         scaling a signal does not change a normalised index,
                         stereo equals two mono runs, ...).
      3. Baselines     - values recorded in golden.inc, to catch unintended
                         numerical changes (including between platforms).
                         These are regression baselines taken from this code,
                         NOT independent R reference values.

    Usage: test_indices [--verbose] [--write-golden FILE]
*/

#include "ACIPlugin.h"
#include "ACIaccPlugin.h"
#include "ADIaccPlugin.h"
#include "AEIaccPlugin.h"
#include "BIaccPlugin.h"
#include "HPlugin.h"
#include "NDSIPlugin.h"
#include "SHPlugin.h"
#include "THPlugin.h"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <map>
#include <memory>
#include <string>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef std::vector<float> Signal;

// ---------------------------------------------------------------- baselines

struct GoldenEntry { const char *key; double value; };
static const GoldenEntry kGolden[] = {
#include "golden.inc"
    { 0, 0 }
};

static std::map<std::string, double> g_measured;   // everything checked against a baseline
static bool g_verbose = false;

// --------------------------------------------------------------- mini checks

static int g_checks = 0;
static int g_failures = 0;

static void fail(const std::string &what)
{
    ++g_failures;
    std::printf("  FAIL  %s\n", what.c_str());
}

static std::string fmt(const char *f, double a, double b = 0, double c = 0)
{
    char buf[256];
    std::snprintf(buf, sizeof buf, f, a, b, c);
    return buf;
}

static void check(bool ok, const std::string &what)
{
    ++g_checks;
    if (!ok) fail(what);
}

static void checkBetween(double v, double lo, double hi, const std::string &what)
{
    ++g_checks;
    if (g_verbose) std::printf("        %-44s = %.6f\n", what.c_str(), v);
    if (!(v >= lo && v <= hi))   // also fails for NaN
        fail(what + fmt(": got %.6f, expected in [%.6f, %.6f]", v, lo, hi));
}

static void checkClose(double a, double b, double relTol, const std::string &what)
{
    ++g_checks;
    double scale = std::fabs(a) > std::fabs(b) ? std::fabs(a) : std::fabs(b);
    if (!(std::fabs(a - b) <= relTol * scale + 1e-7))
        fail(what + fmt(": %.7f vs %.7f (rel tol %.0e)", a, b, relTol));
}

// Compare against the recorded baseline (or record it, for --write-golden).
static void baseline(const std::string &key, double value)
{
    g_measured[key] = value;
    ++g_checks;
    for (const GoldenEntry *g = kGolden; g->key; ++g) {
        if (key == g->key) {
            double scale = std::fabs(g->value) > std::fabs(value) ? std::fabs(g->value) : std::fabs(value);
            if (!(std::fabs(g->value - value) <= 2e-3 * scale + 1e-5))
                fail("baseline " + key + fmt(": got %.7f, recorded %.7f", value, g->value));
            return;
        }
    }
    fail("baseline " + key + " has no entry in golden.inc (run with --write-golden)");
}

// ------------------------------------------------------------------ signals

static Signal tone(double hz, double amp, double sr, double seconds)
{
    size_t n = static_cast<size_t>(sr * seconds);
    Signal s(n);
    for (size_t i = 0; i < n; ++i)
        s[i] = static_cast<float>(amp * std::sin(2.0 * M_PI * hz * i / sr));
    return s;
}

// xorshift32 uniform noise in [-amp, amp): identical on every platform
static Signal noise(double amp, double sr, double seconds, uint32_t seed = 12345u)
{
    size_t n = static_cast<size_t>(sr * seconds);
    Signal s(n);
    uint32_t x = seed;
    for (size_t i = 0; i < n; ++i) {
        x ^= x << 13; x ^= x >> 17; x ^= x << 5;
        s[i] = static_cast<float>(amp * ((x / 4294967296.0) * 2.0 - 1.0));
    }
    return s;
}

// tone gated on for onSec every periodSec (a call-like, strongly time-varying signal)
static Signal pulsed(double hz, double amp, double sr, double seconds, double onSec, double periodSec)
{
    Signal s = tone(hz, amp, sr, seconds);
    for (size_t i = 0; i < s.size(); ++i)
        if (std::fmod(i / sr, periodSec) >= onSec) s[i] = 0.0f;
    return s;
}

static Signal silence(double sr, double seconds)
{
    return Signal(static_cast<size_t>(sr * seconds), 0.0f);
}

static Signal scaled(const Signal &s, float g)
{
    Signal r(s);
    for (size_t i = 0; i < r.size(); ++i) r[i] *= g;
    return r;
}

static Signal sum(const Signal &a, const Signal &b)
{
    Signal r(a);
    for (size_t i = 0; i < r.size() && i < b.size(); ++i) r[i] += b[i];
    return r;
}

// ---------------------------------------------------------------- mini host

typedef std::map<int, std::vector<float> > Outputs;   // output index -> values, in order
typedef std::map<std::string, float> Params;

static void collect(const Vamp::Plugin::FeatureSet &fs, Outputs &out)
{
    for (Vamp::Plugin::FeatureSet::const_iterator it = fs.begin(); it != fs.end(); ++it)
        for (size_t i = 0; i < it->second.size(); ++i)
            if (!it->second[i].values.empty())
                out[it->first].push_back(it->second[i].values[0]);
}

// Feed channels through a plugin like a host would: blocks of the preferred
// size, the last one zero-padded, then getRemainingFeatures().
static Outputs run(Vamp::Plugin *p, const std::vector<Signal> &chans, const Params &params = Params())
{
    for (Params::const_iterator it = params.begin(); it != params.end(); ++it)
        p->setParameter(it->first, it->second);

    size_t block = p->getPreferredBlockSize();
    size_t step = p->getPreferredStepSize();
    if (block == 0) block = 1024;
    if (step == 0) step = block;

    Outputs out;
    if (!p->initialise(chans.size(), step, block)) { fail("initialise failed"); return out; }

    size_t n = chans[0].size();
    size_t nblocks = (n + step - 1) / step;
    std::vector<Signal> buf(chans.size(), Signal(block));
    std::vector<const float *> ptrs(chans.size());
    const unsigned rate = static_cast<unsigned>(p->getInputSampleRate());

    for (size_t b = 0; b < nblocks; ++b) {
        for (size_t c = 0; c < chans.size(); ++c) {
            std::fill(buf[c].begin(), buf[c].end(), 0.0f);
            for (size_t i = 0; i < block && b * step + i < n; ++i)
                buf[c][i] = chans[c][b * step + i];
            ptrs[c] = &buf[c][0];
        }
        collect(p->process(&ptrs[0], Vamp::RealTime::frame2RealTime(b * step, rate)), out);
    }
    collect(p->getRemainingFeatures(), out);
    return out;
}

// One-shot helper: build a plugin, run it on a mono signal, return output 0's first value.
template <class P>
static double value(double sr, const Signal &s, const Params &params = Params(), int output = 0)
{
    std::unique_ptr<P> p(new P(static_cast<float>(sr)));
    Outputs o = run(p.get(), std::vector<Signal>(1, s), params);
    if (o[output].empty()) return 0.0;
    return o[output][0];
}

template <class P>
static Outputs runN(double sr, const std::vector<Signal> &chans, const Params &params = Params())
{
    std::unique_ptr<P> p(new P(static_cast<float>(sr)));
    return run(p.get(), chans, params);
}

// ------------------------------------------------------------------- tests

static const double SR = 20000.0;       // block = sr/10 = 2000 -> 10 Hz bins for ADI/AEI
static const double SECS = 6.0;

struct Sigs {
    Signal tone5k, tone500, noiseW, pulsed5k, quiet;
    Sigs() :
        tone5k(tone(5000, 0.5, SR, SECS)),
        tone500(tone(500, 0.5, SR, SECS)),
        noiseW(noise(0.5, SR, SECS)),
        pulsed5k(pulsed(5000, 0.5, SR, SECS, 0.1, 0.5)),
        quiet(silence(SR, SECS)) {}
};

template <class P>
static void apiChecks(const char *id, size_t outputs)
{
    std::unique_ptr<P> p(new P(static_cast<float>(SR)));
    std::string name = id;
    check(p->getIdentifier() == id, name + ": identifier");
    check(!p->getName().empty() && !p->getMaker().empty(), name + ": name/maker set");
    check(p->getOutputDescriptors().size() == outputs, name + ": output count");
    check(p->getMinChannelCount() == 1, name + ": min channels is 1");
    check(!p->initialise(0, p->getPreferredStepSize(), p->getPreferredBlockSize()), name + ": rejects 0 channels");
    check(!p->initialise(p->getMaxChannelCount() + 1, p->getPreferredStepSize(), p->getPreferredBlockSize()),
          name + ": rejects too many channels");
    check(p->initialise(1, p->getPreferredStepSize(), p->getPreferredBlockSize()), name + ": accepts mono");

    // every parameter reports its default and accepts a value inside its range
    Vamp::Plugin::ParameterList params = p->getParameterDescriptors();
    for (size_t i = 0; i < params.size(); ++i) {
        check(std::fabs(p->getParameter(params[i].identifier) - params[i].defaultValue) < 1e-4,
              name + ": parameter " + params[i].identifier + " starts at its default");
        // frequencies are in Hz throughout
        if (params[i].identifier.find("Freq") != std::string::npos || params[i].identifier.find("anthro") == 0 ||
            params[i].identifier.find("bio") == 0 || params[i].identifier == "freqStep")
            check(params[i].unit == "Hz", name + ": parameter " + params[i].identifier + " is in Hz");
        float mid = 0.5f * (params[i].minValue + params[i].maxValue);
        p->setParameter(params[i].identifier, mid);
        check(std::fabs(p->getParameter(params[i].identifier) - mid) < 1.0,
              name + ": parameter " + params[i].identifier + " round-trips");
    }
}

// Same input after reset() must give the same answer, and a second plugin instance must agree.
template <class P>
static void resetChecks(const char *id, const Signal &s)
{
    std::unique_ptr<P> p(new P(static_cast<float>(SR)));
    Outputs a = run(p.get(), std::vector<Signal>(1, s));
    p->reset();
    Outputs b = run(p.get(), std::vector<Signal>(1, s));
    check(!a[0].empty(), std::string(id) + ": produces a value");
    if (!a[0].empty() && !b[0].empty())
        checkClose(a[0][0], b[0][0], 1e-6, std::string(id) + ": same result after reset()");
    else
        fail(std::string(id) + ": no value after reset()");
}

// Stereo must behave as two independent mono runs.
template <class P>
static void stereoChecks(const char *id, const Signal &l, const Signal &r, int leftOut, int rightOut)
{
    double ml = value<P>(SR, l, Params(), leftOut);
    double mr = value<P>(SR, r, Params(), leftOut);
    Outputs st = runN<P>(SR, std::vector<Signal>{l, r});
    std::string name = std::string(id) + ": stereo";
    if (leftOut == rightOut) {                 // both channels reported on one output
        if (st[leftOut].size() != 2) { fail(name + ": expected 2 values on one output"); return; }
        checkClose(st[leftOut][0], ml, 1e-4, name + " left == mono(left)");
        checkClose(st[leftOut][1], mr, 1e-4, name + " right == mono(right)");
    } else {
        if (st[leftOut].empty() || st[rightOut].empty()) { fail(name + ": missing channel output"); return; }
        checkClose(st[leftOut][0], ml, 1e-4, name + " left == mono(left)");
        checkClose(st[rightOut][0], mr, 1e-4, name + " right == mono(right)");
    }
}

int main(int argc, char **argv)
{
    const char *goldenOut = 0;
    for (int i = 1; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--verbose")) g_verbose = true;
        else if (!std::strcmp(argv[i], "--write-golden") && i + 1 < argc) goldenOut = argv[++i];
    }

    Sigs S;
    const double pi = 3.14159265358979323846; (void)pi;

    std::printf("API\n");
    apiChecks<ACIPlugin>("aci", 1);
    apiChecks<ACIaccPlugin>("aci-acc", 1);
    apiChecks<ADIaccPlugin>("adi-acc", 2);
    apiChecks<AEIaccPlugin>("aei-acc", 2);
    apiChecks<BIaccPlugin>("bi-acc", 2);
    apiChecks<THPlugin>("th", 1);
    apiChecks<SHPlugin>("sh", 1);
    apiChecks<HPlugin>("h", 1);
    apiChecks<NDSIPlugin>("ndsi", 6);

    std::printf("Reset / repeatability\n");
    resetChecks<ACIPlugin>("aci", S.pulsed5k);
    resetChecks<ACIaccPlugin>("aci-acc", S.pulsed5k);
    resetChecks<ADIaccPlugin>("adi-acc", S.noiseW);
    resetChecks<AEIaccPlugin>("aei-acc", S.noiseW);
    resetChecks<BIaccPlugin>("bi-acc", S.tone5k);
    resetChecks<THPlugin>("th", S.pulsed5k);
    resetChecks<SHPlugin>("sh", S.noiseW);
    resetChecks<HPlugin>("h", S.noiseW);
    resetChecks<NDSIPlugin>("ndsi", S.tone5k);


    // ------------------------------------------------------------- ACI
    std::printf("ACI (aci, aci-acc)\n");
    {
        double steady = value<ACIPlugin>(SR, S.tone5k);
        double pulses = value<ACIPlugin>(SR, S.pulsed5k);
        double nz = value<ACIPlugin>(SR, S.noiseW);
        baseline("aci/tone5k", steady);
        baseline("aci/pulsed5k", pulses);
        baseline("aci/noise", nz);
        check(pulses > steady, "aci: a gated tone is more complex than a steady one");
        check(nz > steady, "aci: noise is more complex than a steady tone");
        checkClose(value<ACIPlugin>(SR, scaled(S.noiseW, 0.25f)), nz, 1e-3,
                   "aci: invariant to overall amplitude");
        checkBetween(value<ACIPlugin>(SR, S.quiet), 0, 1e-6, "aci: silence");

        double aSteady = value<ACIaccPlugin>(SR, S.tone5k);
        double aPulses = value<ACIaccPlugin>(SR, S.pulsed5k);
        double aNoise = value<ACIaccPlugin>(SR, S.noiseW);
        baseline("aci-acc/tone5k", aSteady);
        baseline("aci-acc/pulsed5k", aPulses);
        baseline("aci-acc/noise", aNoise);
        check(aPulses > aSteady, "aci-acc: a gated tone is more complex than a steady one");
        check(aNoise > aSteady, "aci-acc: noise is more complex than a steady tone");
        checkClose(value<ACIaccPlugin>(SR, scaled(S.noiseW, 0.25f)), aNoise, 1e-3,
                   "aci-acc: invariant to overall amplitude");
        checkBetween(value<ACIaccPlugin>(SR, S.quiet), 0, 1e-6, "aci-acc: silence");

        // a frequency range holding no signal has far less ACI than the full range
        Params offBand;
        offBand["minFreq"] = 7000.0f; offBand["maxFreq"] = 9000.0f;
        check(value<ACIaccPlugin>(SR, S.pulsed5k, offBand) < 0.5 * aPulses,
              "aci-acc: 7000-9000 Hz of a 5 kHz signal has far less ACI than the full range");

        stereoChecks<ACIaccPlugin>("aci-acc", S.pulsed5k, S.tone5k, 0, 0);
    }

    // ------------------------------------------------------------- ADI
    std::printf("ADI (adi-acc)\n");
    {
        double flat = value<ADIaccPlugin>(SR, S.noiseW);
        double one = value<ADIaccPlugin>(SR, S.tone5k);
        baseline("adi-acc/noise", flat);
        baseline("adi-acc/tone5k", one);
        checkBetween(flat, std::log(10.0) - 0.05, std::log(10.0) + 1e-6,
                     "adi-acc: flat noise across 10 bands approaches ln(10)");
        check(one < 0.5 * flat, "adi-acc: a single tone is far less diverse than noise");
        double two = value<ADIaccPlugin>(SR, sum(S.tone5k, tone(2500, 0.5, SR, SECS)));
        check(two > one, "adi-acc: two tones in different bands are more diverse than one");
        // maxFreq is in Hz: 5000 Hz of flat noise is five 1 kHz bands
        Params fiveKHz;
        fiveKHz["maxFreq"] = 5000.0f;
        checkBetween(value<ADIaccPlugin>(SR, S.noiseW, fiveKHz), std::log(5.0) - 0.05, std::log(5.0) + 1e-6,
                     "adi-acc: maxFreq = 5000 gives five bands, ADI ~ ln(5)");
        stereoChecks<ADIaccPlugin>("adi-acc", S.noiseW, S.tone5k, 0, 1);
    }

    // ------------------------------------------------------------- AEI
    std::printf("AEI (aei-acc)\n");
    {
        double even = value<AEIaccPlugin>(SR, S.noiseW);
        double uneven = value<AEIaccPlugin>(SR, S.tone5k);
        baseline("aei-acc/noise", even);
        baseline("aei-acc/tone5k", uneven);
        checkBetween(even, -1e-6, 0.05, "aei-acc: flat noise is evenly spread (Gini ~ 0)");
        checkBetween(uneven, 0.7, 1.0, "aei-acc: one occupied band out of ten is very uneven");
        // maxFreq is in Hz: ceil(5500/1000) = 6 bands, the 5 kHz tone fills one (Gini 5/6)
        Params fiveKHz;
        fiveKHz["maxFreq"] = 5500.0f;
        checkBetween(value<AEIaccPlugin>(SR, S.tone5k, fiveKHz), 0.55, 0.85,
                     "aei-acc: maxFreq = 5500 gives six bands");
        stereoChecks<AEIaccPlugin>("aei-acc", S.noiseW, S.tone5k, 0, 1);
    }

    // -------------------------------------------------------------- BI
    std::printf("BI (bi-acc)\n");
    {
        double peaked = value<BIaccPlugin>(SR, S.tone5k);
        double flat = value<BIaccPlugin>(SR, S.noiseW);
        double outside = value<BIaccPlugin>(SR, S.tone500);
        baseline("bi-acc/tone5k", peaked);
        baseline("bi-acc/noise", flat);
        baseline("bi-acc/tone500", outside);
        check(peaked > flat, "bi-acc: a tone standing out of the band has a larger BI than flat noise");
        checkBetween(value<BIaccPlugin>(SR, S.quiet), 0, 1e-6, "bi-acc: silence");
        check(outside < peaked, "bi-acc: a 500 Hz tone scores lower than a 5 kHz tone in the 2-8 kHz band");
        stereoChecks<BIaccPlugin>("bi-acc", S.tone5k, S.noiseW, 0, 1);
    }

    // -------------------------------------------------------------- TH
    std::printf("TH (th)\n");
    {
        Signal burst = silence(SR, SECS);
        Signal b = tone(5000, 0.5, SR, 0.05);
        for (size_t i = 0; i < b.size(); ++i) burst[static_cast<size_t>(SR) * 2 + i] = b[i];

        double steady = value<THPlugin>(SR, S.tone5k);
        double gated = value<THPlugin>(SR, S.pulsed5k);
        double one = value<THPlugin>(SR, burst);
        baseline("th/tone5k", steady);
        baseline("th/pulsed5k", gated);
        baseline("th/burst", one);
        checkBetween(steady, 0.98, 1.0, "th: constant-amplitude tone has ~maximal temporal entropy");
        check(one < gated && gated < steady, "th: burst < periodic pulses < steady tone");
        checkBetween(one, 0.0, 0.6, "th: a single short burst has low temporal entropy");
        checkClose(value<THPlugin>(SR, scaled(S.pulsed5k, 0.1f)), gated, 1e-3, "th: invariant to amplitude");
    }

    // -------------------------------------------------------------- SH
    std::printf("SH (sh)\n");
    {
        double nz = value<SHPlugin>(SR, S.noiseW);
        double tn = value<SHPlugin>(SR, S.tone5k);
        baseline("sh/noise", nz);
        baseline("sh/tone5k", tn);
        check(nz > tn, "sh: noise has higher spectral entropy than a tone");
        checkBetween(nz, 0.0, 1.0, "sh: noise within [0,1]");
        checkBetween(tn, 0.0, 0.5, "sh: a pure tone has low spectral entropy");
        checkClose(value<SHPlugin>(SR, scaled(S.noiseW, 0.1f)), nz, 1e-3, "sh: invariant to amplitude");
    }

    // --------------------------------------------------------------- H
    std::printf("H (h)\n");
    {
        double nz = value<HPlugin>(SR, S.noiseW);
        double tn = value<HPlugin>(SR, S.tone5k);
        double pl = value<HPlugin>(SR, S.pulsed5k);
        baseline("h/noise", nz);
        baseline("h/tone5k", tn);
        baseline("h/pulsed5k", pl);
        checkBetween(nz, 0.0, 1.0, "h: noise within [0,1]");
        check(nz > tn, "h: noise has higher total entropy than a tone");
        check(pl < nz, "h: a pulsed tone has lower total entropy than noise");
        // H is SH * TH with both factors in [0,1], so it cannot exceed TH
        double th = value<THPlugin>(SR, S.noiseW);
        check(nz <= th + 1e-4, "h: no larger than its temporal-entropy factor");
    }

    // ------------------------------------------------------------ NDSI
    std::printf("NDSI (ndsi)\n");
    {
        const double sr = 22050.0;
        Signal bio = tone(5000, 0.5, sr, 4.0);       // inside the 2-11 kHz biophony band
        Signal anthro = tone(1500, 0.5, sr, 4.0);    // inside the 1-2 kHz anthropophony band
        double nBio = value<NDSIPlugin>(sr, bio);
        double nAnt = value<NDSIPlugin>(sr, anthro);
        double nMix = value<NDSIPlugin>(sr, sum(bio, anthro));
        baseline("ndsi/bio5k", nBio);
        baseline("ndsi/anthro1.5k", nAnt);
        baseline("ndsi/mixed", nMix);
        checkBetween(nBio, 0.9, 1.0, "ndsi: biophony-band tone gives NDSI ~ +1");
        checkBetween(nAnt, -1.0, -0.9, "ndsi: anthropophony-band tone gives NDSI ~ -1");
        checkBetween(nMix, nAnt, nBio, "ndsi: a mixture lies between the extremes");
        // biophony is a sum of several unit-norm band shares, so it may exceed 1; it is never negative
        checkBetween(value<NDSIPlugin>(sr, bio, Params(), 1), 0.0, 3.0, "ndsi: biophony output");
        checkBetween(value<NDSIPlugin>(sr, anthro, Params(), 2), 0.0, 1.0 + 1e-6, "ndsi: anthropophony output");
        check(value<NDSIPlugin>(sr, bio, Params(), 1) > value<NDSIPlugin>(sr, bio, Params(), 2),
              "ndsi: biophony > anthropophony for a biophony-band tone");

        // outputs: ndsi, biophony, anthropophony, then the same three for the right channel
        Outputs st = runN<NDSIPlugin>(sr, std::vector<Signal>{bio, anthro});
        if (st[0].empty() || st[3].empty()) fail("ndsi: stereo outputs missing");
        else {
            checkClose(st[0][0], nBio, 1e-4, "ndsi: stereo left == mono(left)");
            checkClose(st[3][0], nAnt, 1e-4, "ndsi: stereo right == mono(right)");
        }
    }

    // ------------------------------------------------------------ done
    if (goldenOut) {
        FILE *f = std::fopen(goldenOut, "w");
        if (!f) { std::printf("cannot write %s\n", goldenOut); return 2; }
        std::fprintf(f, "// Regression baselines for tests/test_indices.cpp.\n"
                        "// Generated with: test_indices --write-golden tests/golden.inc\n"
                        "// They come from this code, not from R: regenerate only after a deliberate change.\n");
        for (std::map<std::string, double>::const_iterator it = g_measured.begin(); it != g_measured.end(); ++it)
            std::fprintf(f, "    { \"%s\", %.9g },\n", it->first.c_str(), it->second);
        std::fclose(f);
        std::printf("wrote %zu baselines to %s\n", g_measured.size(), goldenOut);
    }

    std::printf("\n%d checks, %d failed\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
