/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    Vamp

    An API for audio analysis and feature extraction plugins.

    Centre for Digital Music, Queen Mary, University of London.
    Copyright 2006-2012 Chris Cannam and QMUL.
  
    Permission is hereby granted, free of charge, to any person
    obtaining a copy of this software and associated documentation
    files (the "Software"), to deal in the Software without
    restriction, including without limitation the rights to use, copy,
    modify, merge, publish, distribute, sublicense, and/or sell copies
    of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be
    included in all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR
    ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
    CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
    WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

    Except as contained in this notice, the names of the Centre for
    Digital Music; Queen Mary, University of London; and Chris Cannam
    shall not be used in advertising or otherwise to promote the sale,
    use or other dealings in this Software without prior written
    authorization.
*/

#include <vamp-sdk/FFT.h>

/* Ensure pocketfft header does not enable multithreading support
    (prevents pulling in <thread> / pthread dependencies such as libwinpthread). */
#ifndef POCKETFFT_NO_MULTITHREADING
#define POCKETFFT_NO_MULTITHREADING 1
#endif
#include "pocketfft_hdronly.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <vector>
#include <complex>

#if ( VAMP_SDK_MAJOR_VERSION != 2 || VAMP_SDK_MINOR_VERSION != 10 )
#error Unexpected version of Vamp SDK header included
#endif

_VAMP_SDK_PLUGSPACE_BEGIN(FFT.cpp)

namespace Vamp {

void
FFT::forward(unsigned int un,
	     const double *ri, const double *ii,
	     double *ro, double *io)
{
    size_t n = un;
    std::vector<std::complex<double>> data(n);
    for (size_t i = 0; i < n; ++i) {
        data[i] = std::complex<double>(ri[i], ii ? ii[i] : 0.0);
    }

    pocketfft::shape_t shape = {n};
    pocketfft::stride_t stride = {sizeof(std::complex<double>)};
    pocketfft::shape_t axes = {0};
    
    pocketfft::c2c(shape, stride, stride, axes, true, data.data(), data.data(), 1.0);

    for (size_t i = 0; i < n; ++i) {
        ro[i] = data[i].real();
        io[i] = data[i].imag();
    }
}

void
FFT::inverse(unsigned int un,
	     const double *ri, const double *ii,
	     double *ro, double *io)
{
    size_t n = un;
    std::vector<std::complex<double>> data(n);
    for (size_t i = 0; i < n; ++i) {
        data[i] = std::complex<double>(ri[i], ii ? ii[i] : 0.0);
    }

    pocketfft::shape_t shape = {n};
    pocketfft::stride_t stride = {sizeof(std::complex<double>)};
    pocketfft::shape_t axes = {0};
    
    pocketfft::c2c(shape, stride, stride, axes, false, data.data(), data.data(), 1.0);

    double scale = 1.0 / double(n);
    for (size_t i = 0; i < n; ++i) {
        ro[i] = data[i].real() * scale;
        io[i] = data[i].imag() * scale;
    }
}

class FFTComplex::D
{
public:
    D(int n) : m_n(n) { }
    ~D() { }

    void forward(const double *ci, double *co) {
        size_t n = m_n;
        std::vector<std::complex<double>> data(n);
        for (size_t i = 0; i < n; ++i) {
            data[i] = std::complex<double>(ci[i*2], ci[i*2+1]);
        }
        
        pocketfft::shape_t shape = {n};
        pocketfft::stride_t stride = {sizeof(std::complex<double>)};
        pocketfft::shape_t axes = {0};
        
        pocketfft::c2c(shape, stride, stride, axes, true, data.data(), data.data(), 1.0);
        
        for (size_t i = 0; i < n; ++i) {
            co[i*2] = data[i].real();
            co[i*2+1] = data[i].imag();
        }
    }

    void inverse(const double *ci, double *co) {
        size_t n = m_n;
        std::vector<std::complex<double>> data(n);
        for (size_t i = 0; i < n; ++i) {
            data[i] = std::complex<double>(ci[i*2], ci[i*2+1]);
        }
        
        pocketfft::shape_t shape = {n};
        pocketfft::stride_t stride = {sizeof(std::complex<double>)};
        pocketfft::shape_t axes = {0};
        
        pocketfft::c2c(shape, stride, stride, axes, false, data.data(), data.data(), 1.0);
        
        double scale = 1.0 / double(n);
        for (size_t i = 0; i < n; ++i) {
            co[i*2] = data[i].real() * scale;
            co[i*2+1] = data[i].imag() * scale;
        }
    }
    
private:
    size_t m_n;
};

FFTComplex::FFTComplex(unsigned int n) :
    m_d(new D(n))
{
}

FFTComplex::~FFTComplex()
{
    delete m_d;
}

void
FFTComplex::forward(const double *ci, double *co)
{
    m_d->forward(ci, co);
}

void
FFTComplex::inverse(const double *ci, double *co)
{
    m_d->inverse(ci, co);
}

class FFTReal::D
{
public:
    D(int n) : m_n(n) { }
    ~D() { }

    void forward(const double *ri, double *co) {
        size_t n = m_n;
        std::vector<double> in(n);
        for(size_t i=0; i<n; ++i) in[i] = ri[i];
        
        std::vector<std::complex<double>> out(n/2 + 1);
        
        pocketfft::shape_t shape_in = {n};
        pocketfft::stride_t stride_in = {sizeof(double)};
        pocketfft::stride_t stride_out = {sizeof(std::complex<double>)};
        size_t axis = 0;
        
        pocketfft::r2c(shape_in, stride_in, stride_out, axis, true, in.data(), out.data(), 1.0);
        
        for (size_t i = 0; i < n/2 + 1; ++i) {
            co[i*2] = out[i].real();
            co[i*2+1] = out[i].imag();
        }
    }

    void inverse(const double *ci, double *ro) {
        size_t n = m_n;
        std::vector<std::complex<double>> in(n/2 + 1);
        
        for (size_t i = 0; i < n/2 + 1; ++i) {
            in[i] = std::complex<double>(ci[i*2], ci[i*2+1]);
        }
        
        std::vector<double> out(n);
        
        pocketfft::shape_t shape_out = {n};
        pocketfft::stride_t stride_in = {sizeof(std::complex<double>)};
        pocketfft::stride_t stride_out = {sizeof(double)};
        size_t axis = 0;
        
        pocketfft::c2r(shape_out, stride_in, stride_out, axis, false, in.data(), out.data(), 1.0);
        
        double scale = 1.0 / double(n);
        for(size_t i=0; i<n; ++i) {
            ro[i] = out[i] * scale;
        }
    }
    
private:
    size_t m_n;
};

FFTReal::FFTReal(unsigned int n) :
    m_d(new D(n))
{
}

FFTReal::~FFTReal()
{
    delete m_d;
}

void
FFTReal::forward(const double *ri, double *co)
{
    m_d->forward(ri, co);
}

void
FFTReal::inverse(const double *ci, double *ro)
{
    m_d->inverse(ci, ro);
}

}

_VAMP_SDK_PLUGSPACE_END(FFT.cpp)
