// Checks that every plugin reports its declared default for every parameter
// straight after construction, before setParameter is ever called.
#include <cstdio>
#include <memory>
#include "ACIPlugin.h"
#include "ACIaccPlugin.h"
#include "ADIaccPlugin.h"
#include "AEIaccPlugin.h"
#include "BIaccPlugin.h"
#include "THPlugin.h"
#include "SHPlugin.h"
#include "HPlugin.h"
#include "NDSIPlugin.h"

static int failures = 0;

static void check(Vamp::Plugin *p)
{
    Vamp::Plugin::ParameterList params = p->getParameterDescriptors();
    for (size_t i = 0; i < params.size(); ++i) {
        const Vamp::Plugin::ParameterDescriptor &d = params[i];
        float got = p->getParameter(d.identifier);
        bool ok = got == d.defaultValue;
        if (!ok) ++failures;
        printf("%-4s %-12s %-12s default %-10g get %-10g\n", ok ? "ok" : "FAIL",
               p->getIdentifier().c_str(), d.identifier.c_str(), d.defaultValue, got);
    }
}

template <class P> static void run(float rate)
{
    std::unique_ptr<P> p(new P(rate));
    check(p.get());
}

int main()
{
    const float rates[] = { 44100.f, 48000.f, 96000.f };
    for (float r : rates) {
        printf("-- sample rate %g\n", r);
        run<ACIPlugin>(r); run<ACIaccPlugin>(r); run<ADIaccPlugin>(r);
        run<AEIaccPlugin>(r); run<BIaccPlugin>(r); run<THPlugin>(r);
        run<SHPlugin>(r); run<HPlugin>(r); run<NDSIPlugin>(r);
    }
    printf("%d failure(s)\n", failures);
    return failures ? 1 : 0;
}
