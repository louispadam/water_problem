#include "fftw3.h"                   // fftw3 headers
#include "../DoubleVectorNd/SCC_DoubleVector1d.h"    // Array Class Header


#ifndef SCC_FFTW3_COS_1D_
#define SCC_FFTW3_COS_1D_


//
//####################################################################
//                fftw3_cos1d.h : Version Nov 11, 2015
//####################################################################
/**<p>

    This class provides forward and inverse cos transforms of one dimensional
    data stored in GridFunction1d instances or data stored in
    DoubleVector1d instances. Normalization and scaling is performed so that the
    transforms computed are close analogs to the mathematical sin transforms of
    functions defined over intervals of size LX.

    The transforms are implemented using the r2r_1d(...) routine from the
    FFTW3 Fast Fourier Transform library (http://www.fftw.org/) with the
    transform parameter FFTW_REDFT00. The FFTW3 headers and library files
    must be available to create any executable that uses this class.

    This transform interface is designed be used with the data values of
    discrete functions whose values are associated with a uniform grid
    discretization consisting of nx grid panels spanning a
    domain of size LX.

    When DoubleVector1d's are used for input and/or output their size
    is (nx+1); the number of data values associated with
    such discretizations.

    For consistency with other FFTW_InterfaceNd classes, in the
    construction or initialization of class instances, the number
    of panels nx of a discretization is specified.

    For computational efficiency nx should be chosen to be a product of small
    primes (if possible).


    What's computed assuming the input data is specified as a
    GridFunction1d instance initialized with the parameters
    (nx, xMin, xMax):

    ===============================================================

    If k is a wave number then GridFunction1d instance with values given by

    d_i  = sqrt(2.0/LX)*cos(k*pi*(x_i-xMin))/LX)   (when k is in the range 1 <= k <= nx-1)

    d_i  = sqrt(1.0/LX)*cos(k*pi*(x_i-xMin))/LX)   (when k = 0  and k = nx)

    where x_i = xMin + i*hx where i =  0,...,nx and hx = LX/nx

    gets mapped by the forward cosine transform to DoubleVector1d instance of
    size (nx+1), f_hat, with 

    fhat = (0,0,...,1,..,0)

    where the non-zero value entry occurs in the kth position.

    The inverse transform maps the DoubleVector1d instance with

    f_hat = (0,0, ...,1,0,...,0)

    with a value of 1 in the kth location 1 <= k <= nx-1 to  
    containing the values

    d  = sqrt(2.0/LX)*cos(k*pi*(x-xMin))/LX)

    evaluated at the points x = xMin + i*hx where i =  0,...,nx and hx = LX/nx.

    When fhat has a value of 1 in the 0th or nx'th position it is mapped
    to a GridFunction1d instance containing the values

    d  = sqrt(1.0/LX)*cos(k*pi*(x_i-xMin))/LX)  with x = i/nx and i = 0,...,nx

    evaluated at the point x = xMin + i*hx where i =  0,...,nx and hx = LX/nx.

    After computing the transform of the data, one often wants to
    work with the discrete Fourier coefficients. Typically this is done by
    looping over the wave numbers, and then obtaining the coefficient by
    computing the correct offset. For example, one uses code constructs
    as indicated by the fragment below.

    //
    // Loop over cos wave numbers
    //

    for(long kx = 0; kx <= nx; kx++)
    {
        sCoeff  = f_hat(kx);  // (kx)'th cos coefficient is
                *             // the (kx)'th entry of the transform
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

namespace SCC
{
class fftw3_cos1d
{
public:


fftw3_cos1d()
{
    plan     = nullptr;
    in       = nullptr;
    out      = nullptr;
    nx       = 0;
    LX       = 0;
    nSamples = 0;
}

fftw3_cos1d(long nx,double LX = 1.0)
{
    this->nx    = 0;
    in          = nullptr;
    out         = nullptr;
    plan        = nullptr;
    initialize(nx,LX);
}


virtual ~fftw3_cos1d()
{
    if(plan != nullptr)
    {fftw_destroy_plan(plan);}
    
    if(in  != nullptr)  fftw_free(in);
    if(out != nullptr)  fftw_free(out);
}

void initialize()
{
    if(plan != nullptr)
    {fftw_destroy_plan(plan);}
    
    if(in   != nullptr) fftw_free(in);
    if(out  != nullptr) fftw_free(out);

    plan     = nullptr;
    in       = nullptr;
    out      = nullptr;

    nx       = 0;
    LX       = 1.0;
    nSamples = 0;
}

void initialize(long nx,double LX = 1.0)
{
    if((this->nx != nx))
    {
    this->nx       = nx;
    this->nSamples = nx+1;
    if(plan != nullptr)
    {fftw_destroy_plan(plan);}

    if(in  != nullptr) fftw_free(in);
    if(out != nullptr) fftw_free(out);

    in  = (double*) fftw_malloc(sizeof(double) * nSamples);
    out = (double*) fftw_malloc(sizeof(double) * nSamples);
    plan = fftw_plan_r2r_1d(nSamples, in, out, FFTW_REDFT00, FFTW_ESTIMATE);

    if(plan == nullptr)
    {
    throw std::runtime_error("\nXXX Error : required fftw_r2r_1d function not available \nXXX in FFTW library used (likely MKL) ");
    }
    }

    this->LX = LX;
}


// replan() is used to create new plans for multi-threaded instances when
// the number of threads allocated to FFTW3 is changed.

void replan()
{
   if(plan != nullptr) {fftw_destroy_plan(plan);}

   if((in  == nullptr) || (out == nullptr))
   {
   throw std::runtime_error("\nXXX Error : calling replan() before SCC::fft3_sin1d instance initialized.\n ");
   }

   plan = fftw_plan_r2r_1d(nSamples, in, out, FFTW_REDFT00, FFTW_ESTIMATE);

   if(plan == nullptr)
   {
   throw std::runtime_error("\nXXX Error : required fftw_r2r_1d function not available \nXXX in FFTW library used (likely MKL) ");
   }
}

// fftw2d_cos_forward argument sizes:
//
// DoubleVector1d size nx+1
//

void fftw1d_cos_forward(DoubleVector1d&  F)
{
    long k; 
     double scalingfactor;
     
    if(nx != F.getSize()-1)
    { 
    initialize(F.getSize()-1);
    }
    
    //copy input
    
    for(k=0; k < nSamples; k++)
    {
        in[k] = F(k);
    }

    
    fftw_execute(plan);
    

    scalingfactor = sqrt(LX)/(sqrt(2.0)*((double)(nx)));
    
    F(0)          = 0.5*out[0]*sqrt(2.0)*scalingfactor;
    F(nSamples-1) = 0.5*out[nSamples-1]*sqrt(2.0)*scalingfactor;
    for(k=1; k < nSamples-1; k++)
    {
        F(k) = out[k]*scalingfactor;
    }
}

// fftw2d_cos_forward argument sizes:
//
// DoubleVector1d: Initialized to size nx+1
//

void fftw1d_cos_forward(DoubleVector1d&  inF,  DoubleVector1d& outF)
{
    long k; 
     double scalingfactor;

    if(nx != inF.getSize()-1)
    { 
    initialize(inF.getSize()-1);
    }
    
    //copy input
    
    for(k=0; k < nSamples; k++)
    {
        in[k] = inF(k);
    }

    fftw_execute(plan);
    
    scalingfactor = sqrt(LX)/(sqrt(2.0)*((double)(nx)));
    
    outF(0)          = 0.5*out[0]*sqrt(2.0)*scalingfactor;
    outF(nSamples-1) = 0.5*out[nSamples-1]*sqrt(2.0)*scalingfactor;
    for(k=1; k < nSamples-1; k++)
    {
        outF(k) = out[k]*scalingfactor;
    }

}


// fftw2d_cos_inverse argument sizes:
//
// DoubleVector1d: Initialized to size nx+1
//

void fftw1d_cos_inverse(DoubleVector1d&  F)
{
    double scalingfactor;

    if(nx != F.getSize()-1)
    { 
    initialize(F.getSize()-1);
    }

    //copy input and scale
    
    scalingfactor = 1.0/sqrt(2.0*LX);
    for(long k=0; k < nSamples; k++)
    {
        in[k] = F(k)*scalingfactor;
    }

    in[0]          *= sqrt(2.0);
    in[nSamples-1] *= sqrt(2.0);

    fftw_execute(plan);

    for(long k=0; k < nSamples; k++)
    {
    F(k) = out[k];
    }
}

// fftw2d_cos_inverse argument sizes:
//
// DoubleVector1d: Initialized to size nx+1
//

void fftw1d_cos_inverse(DoubleVector1d&  inF,  DoubleVector1d& outF)
{
    double scalingfactor;

    if(nx != inF.getSize()-1)
    { 
    initialize(inF.getSize()-1);
    }

    //copy input
    
    scalingfactor = 1.0/sqrt(2.0*LX);
    for(long k=0; k < nSamples; k++)
    {
        in[k] = inF(k)*scalingfactor;
    }

    in[0]          *= sqrt(2.0);
    in[nSamples-1] *= sqrt(2.0);


    fftw_execute(plan);

    for(long k=0; k < nSamples; k++)
    {
    outF(k) = out[k];
    }
}

private:

    long nx; double LX;
    long nSamples;

    fftw_plan plan;

    double*  in;
    double* out;

};
}
#endif
