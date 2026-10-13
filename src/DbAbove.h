/*
    Ecoacoustic Vamp Plugins

    High-performance implementations of acoustic indices for
    bioacoustics and soundscape ecology analysis.

    (C) Ed Baker 2025. Licensed under GPL (>=2).
*/

// DbAbove.h - "is this magnitude above a threshold in dB?", with log10 only where it can matter.
//
// The plugins count cells whose level db = 20 * log10(magnitude) (db = -200 for magnitudes up to 1e-10) is above a
// threshold in dB. Almost every cell is far from the threshold, so the answer is clear from the magnitude alone:
// the decision is made by comparing it with 10^((threshold +- margin) / 20), and the log10 of the original test is
// evaluated, exactly as before, only for the few magnitudes between those two bounds.
//
// Why the shortcut cannot change an answer: with the margin at 2e-6 dB, a magnitude beyond either bound has a
// true level more than about 1.99e-6 dB from the threshold (the bounds themselves are off by under 1e-10 dB),
// whereas the computed 20 * log10(magnitude) is within a few 1e-14 dB of the true level (log10 is accurate to
// an ulp or two, then one multiplication). So the computed level is on the same side of the threshold.

#ifndef _DB_ABOVE_H_
#define _DB_ABOVE_H_

#include <cmath>

class DbAbove
{
public:
    explicit DbAbove(double thresholdDB) :
        m_threshold(thresholdDB),
        m_bounded(std::isfinite(thresholdDB) && std::fabs(thresholdDB) < 1000.0),
        m_low(0.0), m_high(0.0)
    {
        if (m_bounded) {
            m_high = std::pow(10.0, (thresholdDB + MARGIN_DB) / 20.0);
            m_low = std::pow(10.0, (thresholdDB - MARGIN_DB) / 20.0);
        }
    }

    // Same result as:  double db = -200.0; if (magnitude > 1e-10) db = 20.0 * std::log10(magnitude);  return db > thresholdDB;
    bool operator()(double magnitude) const {
        if (m_bounded && magnitude > 1e-10) {
            if (magnitude > m_high) return true;
            if (magnitude < m_low) return false;
        }
        double db = -200.0;
        if (magnitude > 1e-10) {
            db = 20.0 * std::log10(magnitude);
        }
        return db > m_threshold;
    }

private:
    static constexpr double MARGIN_DB = 2e-6;
    double m_threshold;
    bool m_bounded;
    double m_low;
    double m_high;
};

#endif
