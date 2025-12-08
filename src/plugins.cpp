/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    ACI Vamp Plugin
    
    Implements the Acoustic Complexity Index as described in 
    Pieretti et al. (2011) and mimics seewave::ACI()
*/

#include "vamp/vamp.h"
#include "vamp-sdk/PluginAdapter.h"

#include "ACIPlugin.h"
#include "ACIaccPlugin.h"
#include "ADIaccPlugin.h"
#include "ADIPlugin.h"
#include "AEIPlugin.h"
#include "AEIaccPlugin.h"
#include "BIaccPlugin.h"
#include "BIPlugin.h"
#include "AmplitudeIndexPlugin.h"
#include "THPlugin.h"
#include "SHPlugin.h"
#include "HPlugin.h"
#include "NDSIPlugin.h"

static Vamp::PluginAdapter<ACIPlugin> aciAdapter;
static Vamp::PluginAdapter<ACIaccPlugin> aciAccAdapter;
static Vamp::PluginAdapter<ADIaccPlugin> adiAccAdapter;
static Vamp::PluginAdapter<ADIPlugin> adiAdapter;
static Vamp::PluginAdapter<AEIPlugin> aeiAdapter;
static Vamp::PluginAdapter<AEIaccPlugin> aeiAccAdapter;
static Vamp::PluginAdapter<BIaccPlugin> biAccAdapter;
static Vamp::PluginAdapter<BIPlugin> biAdapter;
static Vamp::PluginAdapter<AmplitudeIndexPlugin> amplitudeIndexAdapter;
static Vamp::PluginAdapter<THPlugin> thAdapter;
static Vamp::PluginAdapter<SHPlugin> shAdapter;
static Vamp::PluginAdapter<HPlugin> hAdapter;
static Vamp::PluginAdapter<NDSIPlugin> ndsiAdapter;

const VampPluginDescriptor *vampGetPluginDescriptor(unsigned int version,
                                                    unsigned int index)
{
    if (version < 1) return 0;

    switch (index) {
    case  0: return aciAdapter.getDescriptor();
    case  1: return amplitudeIndexAdapter.getDescriptor();
    case  2: return thAdapter.getDescriptor();
    case  3: return aciAccAdapter.getDescriptor();
    case  4: return adiAccAdapter.getDescriptor();
    case  5: return adiAdapter.getDescriptor();
    case  6: return aeiAdapter.getDescriptor();
    case  7: return aeiAccAdapter.getDescriptor();
    case  8: return biAccAdapter.getDescriptor();
    case  9: return biAdapter.getDescriptor();
    case 10: return shAdapter.getDescriptor();
    case 11: return hAdapter.getDescriptor();
    case 12: return ndsiAdapter.getDescriptor();
    default: return 0;
    }
}
