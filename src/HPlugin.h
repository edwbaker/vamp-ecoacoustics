#ifndef _H_PLUGIN_H_
#define _H_PLUGIN_H_

#include "SHPlugin.h"
#include "THPlugin.h"

class HPlugin : public SHPlugin
{
public:
    HPlugin(float inputSampleRate);
    virtual ~HPlugin();

    string getIdentifier() const;
    string getName() const;
    string getDescription() const;
    string getMaker() const;
    int getPluginVersion() const;
    string getCopyright() const;

    size_t getPreferredBlockSize() const;
    size_t getPreferredStepSize() const;

    OutputList getOutputDescriptors() const;

    bool initialise(size_t channels, size_t stepSize, size_t blockSize);
    void reset();

    FeatureSet process(const float *const *inputBuffers,
                       Vamp::RealTime timestamp);

    FeatureSet getRemainingFeatures();

protected:
    THPlugin* m_thPlugin;
};

#endif
