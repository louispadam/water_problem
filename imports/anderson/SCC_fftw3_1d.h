//
//####################################################################
//   fftw3_1d.h : Version Feb. 7, 2013
//####################################################################
//
#include "fftw3.h"                                    // fftw3 headers
#include "../DoubleVectorNd/SCC_DoubleVector1d.h"   // Array Class Header
#include "../GridFunctionNd/SCC_GridFunction1d.h"   // GridFunction Class Header

#ifndef SCC_FFTW3_1D_
#define SCC_FFTW3_1D_

namespace SCC
{
//
//####################################################################
//                     fftw3_1d.h : Version Nov. 11, 2015
//####################################################################
/**<p>

    This class provides forward and inverse transforms of one dimensional
    data stored in GridFunction1d instances or data stored in
    DoubleVector1d instances. Normalization and scaling is performed so that the
    transforms computed are close analogs to the mathematical transforms of
    periodic functions defined over intervals of size LX.

    The transforms are implemented using the dft_1d(...) routine from the
    FFTW3 Fast Fourier Transform library (http://www.fftw.org/). The
    FFTW3 headers and library files must be available to create any executable
    that uses this class.

    This transform interface is designed be used with the data values of
    discrete periodic functions associated with a uniform grid discretization
    consisting of nx grid panels spanning a domain of size LX.

    For consistency with other transform interface routines the size of the
    input and output data arrays consists of data of size nx.

    For computational efficiency nx should be chosen to be a product of small
    primes (if possible).

    ===============================================================

    What's computed assuming the input data is specified as a
    GridFunction1d instance initialized with the parameters
    (nx, xMin, xMax):

    If kx is a wave number in the range from -(nx/2) to (nx-1)/2
    (the range calculated with integer division), then two GridFunction1d
    instances whose values are specified by the real and imaginary parts of
    sqrt(1/LX)*exp(2*pi*I*kx*(x-xMin)/LX)

    evaluated at the points x = xMin + i*hx where i = 0,1,...,nx and hx = LX/nx

    gets mapped by the forward transform to the two DoubleVector2d instances
    of size nx, fhat_real and fhat_imag, where 

    fhat_real  = (0,0, ... 1, ... 0)
    fhat_imag  = (0,0, ... 1, ... 0)

    Here the non-zero entry occurs in the kx + (nx/2) location.

    Similarly the inverse transform maps two DoubleVector2d instances

    fhat_real = (0,0, ... 1, ... 0)
    fhat_omag = (0,0, ... 1, ... 0)

    with a value of 1 in the k + (nx/2) location two GridFunction1d instances
    consisting of the real and imaginary parts of

    d = sqrt(1/LX)*exp(2*pi*I*kx*(x-xMin)/LX)

    evaluated at the points  x = xMin + i*hx where i = 0,1,...,nx and hx = LX/nx

    After computing the transform of the data, one often wants to
    work with the discrete Fourier coefficients. Typically this is done by
    looping over the wave numbers with the corresponding coefficient index
    determined using an offset. For example, if one creates the real
    and imginary transform components fhat_real and fhat_imag, then
    the following code fragment is can be used to work with the
    real and imaginary coefficients of this transform:

    Note: (nx/2) and (nx-1)/2 must be computed using integer division!

   //
   // Loop over wave numbers kx
   //

   long kxIndex;

   for(kx = -(nx/2); kx <= (nx-1)/2; k++)
   {
    kxIndex = kx + (nx/2);      // index of kx'th transform coefficient

    rk = fhat_real(kIndex);     // rk =  real coefficient of the kx'th basis function
    ik = fhat_imag(kIndex);     // ik =  imag coefficient of the kx'th basis function
             *
             *
             *

  }
*/
/*
#############################################################################
#
# Copyright 2014-2017 Chris Anderson
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# For a copy of the GNU General Public License see
# <http://www.gnu.org/licenses/>.
#
#############################################################################
*/

class fftw3_1d
{
public:

fftw3_1d()
{
    forwardplan = nullptr;
    inverseplan = nullptr;

    in  = nullptr;
    out = nullptr;
    nx  = 0;
    LX  = 1.0;
}

fftw3_1d(long nx, double LX = 1.0)
{
    this->nx    = 0;
    in          = nullptr;
    out         = nullptr;
    forwardplan = nullptr;
    inverseplan = nullptr;

    initialize(nx,LX);
}

fftw3_1d(const fftw3_1d& DFT)
{
    if(DFT.forwardplan == nullptr)
    {
        in  = nullptr;
        out = nullptr;
        forwardplan = nullptr;
        inverseplan = nullptr;
        initialize();
        return;
    }
     initialize(DFT.nx,DFT.LX);
}

virtual ~fftw3_1d()
{
    if(forwardplan != nullptr)
    { fftw_destroy_plan(forwardplan);}
    
    if(inverseplan != nullptr)
    {fftw_destroy_plan(inverseplan);}

    if(in  != nullptr) {fftw_free(in);}
    if(out != nullptr) {fftw_free(out);}
}

void initialize()
{
    if(forwardplan != nullptr)
    { fftw_destroy_plan(forwardplan);}
    
    if(inverseplan != nullptr)
    {fftw_destroy_plan(inverseplan);}

    
    if(in  != 0) {fftw_free(in);}
    if(out != 0) {fftw_free(out);}

    forwardplan = nullptr;
    inverseplan = nullptr;

    in  = nullptr;
    out = nullptr;
    nx  = 0;
    LX  = 1.0;
}

void initialize(long nx, double LX = 1.0)
{
    if((this->nx != nx))
    {
    this->nx = nx;
    
    if(forwardplan != nullptr)
    { fftw_destroy_plan(forwardplan);}
    
    if(inverseplan != nullptr)
    {fftw_destroy_plan(inverseplan);}

    if(in  != nullptr)  fftw_free(in);
    if(out != nullptr) fftw_free(out);

    in  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx);

    forwardplan = fftw_plan_dft_1d(nx, in, out, FFTW_FORWARD,  FFTW_ESTIMATE);
    inverseplan = fftw_plan_dft_1d(nx, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

    // Note: In the FFTW_BACKWARD plan the "out" argument is transformed 
    //       into the "in" argument when the plan is executed. 
    }

    this->LX = LX;
}


// replan() is used to create new plans for multi-threaded instances when
// the number of threads allocated to FFTW3 is changed.

void replan()
{
   if(forwardplan != nullptr)
   { fftw_destroy_plan(forwardplan);}

   if(inverseplan != nullptr)
   {fftw_destroy_plan(inverseplan);}

   if((in  == nullptr) || (out == nullptr))
   {
   throw std::runtime_error("\nXXX Error : calling replan() before SCC::fft3_1d instance initialized.\n ");
   }

   forwardplan = fftw_plan_dft_1d(nx, in, out, FFTW_FORWARD,  FFTW_ESTIMATE);
   inverseplan = fftw_plan_dft_1d(nx, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
}


// fftw1d_forward argument sizes:
//
// DoubleVector1d: Initialized of size nx

void fftw1d_forward(DoubleVector1d&  inReal,  DoubleVector1d& inImag,
                    DoubleVector1d& outReal,  DoubleVector1d& outImag)
{
    long k; long i;

    if(nx != inReal.getSize())
    { 
    initialize(inReal.getSize());
    }

    //copy input
    
    for(k=0; k < nx; k++)
    {
        in[k][0] = inReal(k);
        in[k][1] = inImag(k);
    }

    fftw_execute(forwardplan);
    
    // extract output and scale

    double scalefactor = sqrt(LX)/nx;

    for(k=0; k < nx; k++)
    {
        i =   (k + (nx/2))%nx;   //NB: uses integer division of odd numbers

        outReal(i) = out[k][0]*scalefactor;
        outImag(i) = out[k][1]*scalefactor;
    }
}

// fftw1d_inverse argument sizes:
//
// DoubleVector1d: Initialized of size nx
//


void fftw1d_inverse (DoubleVector1d&  inReal,  DoubleVector1d& inImag,
                     DoubleVector1d& outReal,  DoubleVector1d& outImag)
{
    if(nx != inReal.getSize())
    {
    initialize(inReal.getSize());
    }

    long i; long k;

    //reorder input

    for(i=0; i < nx; i++)
    {
        k = (i - (nx/2) + nx )%nx ;  //NB: uses integer division of odd numbers
        in[k][0] = inReal(i);
        in[k][1] = inImag(i);
    }

    //transform

    fftw_execute(inverseplan);

    // Capture periodic data

    double scalefactor =  1.0/sqrt(LX);

    for(k=0; k < nx; k++)
    {
        outReal(k) = out[k][0]*scalefactor;
        outImag(k) = out[k][1]*scalefactor;
    }
}


// fftw1d_forward argument sizes:
//
// GridFunction1d: Initialized with nx panels
// DoubleVector1d: Initialized of size nx

void fftw1d_forward(GridFunction1d&  inReal,  GridFunction1d& inImag,
                    DoubleVector1d& outReal,  DoubleVector1d& outImag)
{
    long k; long i;

    // Capture domain size

    this->LX = inReal.getXmax() - inReal.getXmin();

    if(nx != inReal.getXpanelCount())
    {
    initialize(inReal.getXpanelCount(),LX);
    }

    //copy input

    for(k=0; k < nx; k++)
    {
        in[k][0] = inReal(k);
        in[k][1] = inImag(k);
    }

    fftw_execute(forwardplan);

    double scalefactor =  sqrt(LX)/nx;

    for(k=0; k < nx; k++)
    {
        i =   (k + (nx/2))%nx;   //NB: uses integer division of odd numbers

        outReal(i) = out[k][0]*scalefactor;
        outImag(i) = out[k][1]*scalefactor;
    }
}

// fftw1d_inverse argument sizes:
//
// DoubleVector1d: Initialized of size nx
// GridFunction1d: Initialized with nx panels


void fftw1d_inverse(DoubleVector1d&  inReal,  DoubleVector1d& inImag,
                    GridFunction1d& outReal,  GridFunction1d& outImag)
{
    this->LX = outReal.getXmax() - outReal.getXmin();

    if(nx != inReal.getSize())
    { 
    initialize(inReal.getSize(),LX);
    }

    long i; long k; 

    //reorder input
    
    for(i=0; i < nx; i++)
    {
        k = (i - (nx/2) + nx )%nx ;  //NB: uses integer division of odd numbers
        in[k][0] = inReal(i);
        in[k][1] = inImag(i);
    }

    //transform

    fftw_execute(inverseplan);
    
    // Capture periodic data

    double scalefactor =  1.0/sqrt(LX);

    for(k=0; k < nx; k++)
    {
        outReal(k) = out[k][0]*scalefactor;
        outImag(k) = out[k][1]*scalefactor;
    }

    // Enforce periodicity of output GridFunction

    outReal.enforcePeriodicity();
    outImag.enforcePeriodicity();
}

public:

    long nx; double LX;

    fftw_plan forwardplan;
    fftw_plan inverseplan;

    fftw_complex*  in;
    fftw_complex* out;
};
}
#endif
