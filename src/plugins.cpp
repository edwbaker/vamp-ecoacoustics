/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    ACI Vamp Plugin
    
    Implements the Acoustic Complexity Index as described in 
    Pieretti et al. (2011) and mimics seewave::ACI()
*/

#include "vamp/vamp.h"
#include "vamp-sdk/PluginAdapter.h"

#include "ACIPlugin.h"
#include "AmplitudeIndexPlugin.h"
#include "TemporalEntropyPlugin.h"
#include "ACImtPlugin.h"

static Vamp::PluginAdapter<ACIPlugin> aciAdapter;
static Vamp::PluginAdapter<AmplitudeIndexPlugin> amplitudeIndexAdapter;
static Vamp::PluginAdapter<TemporalEntropyPlugin> temporalEntropyAdapter;
static Vamp::PluginAdapter<ACImtPlugin> acimtAdapter;

const VampPluginDescriptor *vampGetPluginDescriptor(unsigned int version,
                                                    unsigned int index)
{
    if (version < 1) return 0;

    switch (index) {
    case  0: return aciAdapter.getDescriptor();
    case  1: return amplitudeIndexAdapter.getDescriptor();
    case  2: return temporalEntropyAdapter.getDescriptor();
    case  3: return acimtAdapter.getDescriptor();
    default: return 0;
    }
}
