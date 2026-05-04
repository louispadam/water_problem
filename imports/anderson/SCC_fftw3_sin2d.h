#include "fftw3.h"
#include "../DoubleVectorNd/SCC_DoubleVector2d.h"    // Array class header
#include "../GridFunctionNd/SCC_GridFunction2d.h"    // Grid function class header

#include <exception>


#ifndef SCC_FFTW3_SIN_2D_
#define SCC_FFTW3_SIN_2D_
//
//####################################################################
//          SCC_fftw3_sin2d.h : Version Nov. 12 2015
//####################################################################
//
/**<p>
    This class provides forward and inverse sin transforms of two dimensional
    data stored in GridFunction2d instances or data stored in
    DoubleVector2d instances. Normalization and scaling is performed so that the
    transforms computed are close analogs to the mathematical sin transforms of
    functions defined over intervals of size LX x LY.

    The transforms are implemented using the r2r_2d(...) routines from the
    FFTW3 Fast Fourier Transform library (http://www.fftw.org/). The
    FFTW3 headers and library files must be available to create any executable
    that uses this class.

    This transform interface is designed be used with the data values of
    discrete functions with homogeneous boundary values associated with
    a uniform grid discretization consisting of nx grid panels by ny grid panels
    spanning a domain of size LX x LY.

    When DoubleVector2d's are used for input and/or output their size
    is (nx-1) x (ny-1); the number of interior data values associated with
    such discretizations.

    For consistency with other FFTW_InterfaceNd classes, in the
    construction or initialization of class instances, the number
    of panels nx and ny of a discretization is specified.

    For computational efficiency nx and ny should be chosen to be a
    product of small primes less than or equal to 13  (if possible).

    ===============================================================

    What's computed assuming the input data is specified as a
    GridFunction2d instance initialized with the parameters
    (nx, xMin, xMax, ny, yMin, yMax):

    If (kx,ky) is a wave number vector with kx, and ky in the ranges

    1 <= kx <=  (nx-1)   and 1 <= ky <= (ny-1)

    then a GridFunction2d instance with values given by

    d  = sqrt(2/LX)*sin(kx*pi*(x_i-xMin)/LX)
        *sqrt(2/LY)*sin(ky*pi*(y_j-yMin)/LY)

    evaluated the points

    x_i = xMin + i*hx, i = 0,1,...,nx-1, hx = LX/nx
    y_j = yMin + j*hy, j = 0,1,...,ny-1, hy = LY/ny

    is mapped by the forward transform to a DoubleVector2d of size
    (nx-1)x*(ny-1), f_hat, where

    f_hat(kx-1, ky-1) = 1

    and all other values are zero.

    Similarly the inverse transform maps the DoubleVector2d of size
    (nx-1)x*(ny-1), f_hat, that is zero value except for a value of 1 in
    (kx-1,ky-1) place to the GridFunction2d instance with values
    corresponding to

    d  = sqrt(2/LX)*sin(kx*pi*(x_i-xMin)/LX)
        *sqrt(2/LY)*sin(ky*pi*(y_j-yMin)/LY)

    evaluated the points

    x = xMin + i*hx, i = 0,1,...,nx-1, hx = LX/nx
    y = yMin + j*hy, j = 0,1,...,ny-1, hy = LY/ny

    After computing the transform of the data, one often wants to
    work with the discrete Fourier coefficients. Typically this is done by
    looping over the wave numbers, and then obtaining the coefficient by
    computing the correct offset. For example, one uses code constructs
    as indicated by the fragment below.

    //
    // Loop over sin wave numbers
    //

    for(long kx = 1; kx <= nx-1; kx++)
    {
    for(long ky = 1; ky <= ny-1; ky++)
    {
        sCoeff  = f_hat(kx-1,ky-1);  // (kx,ky)'th  sin coefficient is
                *                    // the (kx-1,ky-1) entry of the transform
                *
                *
    }}
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
//####################################################################
// Chris Anderson (C) UCLA                               Nov. 12, 2015
//####################################################################
//
namespace SCC
{
class fftw3_sin2d
{
public:

fftw3_sin2d()
{
    plan = nullptr;
    in   = nullptr;
    out  = nullptr;

    nx  = 0;
    ny  = 0;

    LX = 1.0;
    LY = 1.0;

    nSampleX = 0;
    nSampleY = 0;
}

// Initialization values of nx and ny correspond to the number of panels
// of a discretization of a domain of size LX x LY

fftw3_sin2d(long nx, long ny,double LX = 1.0, double LY = 1.0)
{
    this->nx    = 0;
    this->ny    = 0;

    in          = nullptr;
    out         = nullptr;
    plan        = nullptr;

    initialize(nx,ny);
}

virtual ~fftw3_sin2d()
{
    if(plan != nullptr) {fftw_destroy_plan(plan);}

    if(in  != nullptr) fftw_free(in);
    if(out != nullptr) fftw_free(out);

}

void initialize()
{
    if(plan != nullptr)
    {fftw_destroy_plan(plan);};

    if(in  != nullptr) fftw_free(in);
    if(out != nullptr) fftw_free(out);


    plan = nullptr;

    in  = nullptr;
    out = nullptr;

    nx  = 0;
    ny  = 0;

    LX = 1.0;
    LY = 1.0;

    nSampleX = 0;
    nSampleY = 0;
}

void initialize(long nx, long ny,double LX = 1.0, double LY = 1.0)
{
    if((this->nx != nx)||(this->ny != ny))
    {
    this->nx   = nx;
    this->ny   = ny;
    nSampleX   = nx-1;
    nSampleY   = ny-1;
    if(plan != nullptr)
    {fftw_destroy_plan(plan);}

    if(in  != nullptr) fftw_free(in);
    if(out != nullptr) fftw_free(out);

    in  = (double*) fftw_malloc(sizeof(double) * nSampleX*nSampleY);
    out = (double*) fftw_malloc(sizeof(double) * nSampleX*nSampleY);

    plan = fftw_plan_r2r_2d(nSampleX,nSampleY, in, out,FFTW_RODFT00,
           FFTW_RODFT00,FFTW_ESTIMATE);

    if(plan == nullptr)
    {
    throw std::runtime_error("\nXXX Error : required fftw_r2r_2d function not available \nXXX in FFTW library used (likely MKL) ");
    }

    // Efficiency: From the documentation for RODFT00 transform, nSampleX+1 and nSampleY+1
    // should be a product of small primes ==> nx and ny should be product of small primes.
    }

    this->LX = LX;
    this->LY = LY;
}
    

// replan() is used to create new plans for multi-threaded instances when
// the number of threads allocated to FFTW3 is changed.

void replan()
{
   if(plan != nullptr) {fftw_destroy_plan(plan);}

   if((in  == nullptr) || (out == nullptr))
   {
   throw std::runtime_error("\nXXX Error : calling replan() before SCC::fft3_sin2d instance initialized.\n ");
   }

   plan = fftw_plan_r2r_2d(nSampleX,nSampleY, in, out,FFTW_RODFT00, FFTW_RODFT00,FFTW_ESTIMATE);

   if(plan == nullptr)
   {
   throw std::runtime_error("\nXXX Error : required fftw_r2r_2d function not available \nXXX in FFTW library used (likely MKL) ");
   }
}

// fftw2d_sin_forward argument sizes:
//
// DoubleVector2d: Initialized to size nx-1 by ny-1
//
// This operator ignores the perimeter values of the input vector.

void fftw2d_sin_forward(GridFunction2d& inF, DoubleVector2d& outF)
{
    long i,j;

    this->LX = inF.getXmax() - inF.getXmin();
    this->LY = inF.getYmax() - inF.getYmin();

    if((nx != inF.getXpanelCount()) || (ny != inF.getYpanelCount()))
    {
    initialize(inF.getXpanelCount(),inF.getYpanelCount(),LX,LY);
    }

    // Extract input, ignoring perimeter data

    for(i=0; i < nSampleX; i++)
    {
    for(j=0; j < nSampleY; j++)
    {
        in[j+ i*nSampleY] = inF(i+1,j+1);
    }}

    fftw_execute(plan);

    // Capture transform values and scale appropriately

    double scalingfactor = sqrt(LX*LY)/(2.0*((double)(nx)*(double)(ny)));

    for(i=0; i < nSampleX; i++)
    {
    for(j=0; j < nSampleY; j++)
    {
        outF(i,j) = out[j+ i*nSampleY]*scalingfactor;
    }}

}
// fftw2d_sin_inverse argument sizes:
//
// DoubleVector2d: Initialized to size nx-1 by ny-1
//
// This operator ignores the perimeter values of the input vector.

void fftw2d_sin_inverse(DoubleVector2d& inF, GridFunction2d& outF)
{
    long i,j;

    this->LX = outF.getXmax() - outF.getXmin();
    this->LY = outF.getYmax() - outF.getYmin();

    if((nx != inF.getIndex1Size()+1) || (ny != inF.getIndex2Size()+1))
    {
    initialize(inF.getIndex1Size()+1,inF.getIndex2Size()+1,LX,LY);
    }

    for(i=0; i < nSampleX; i++)
    {
    for(j=0; j < nSampleY; j++)
    {
        in[j+ i*nSampleY] = inF(i,j);
    }}

    fftw_execute(plan);

    double scalingfactor = 1.0/(2.0*sqrt(LX*LY));

    for(i=0; i < nSampleX; i++)
    {
    for(j=0; j < nSampleY; j++)
    {
        outF(i+1,j+1) = out[j+ i*nSampleY]*scalingfactor;
    }}

    outF.setBoundaryValues(0.0);
}

// fftw2d_sin_forward argument sizes:
//
// DoubleVector2d size nx-1 by ny-1
//

void fftw2d_sin_forward(DoubleVector2d& inF, DoubleVector2d& outF)
{
    long i,j;

    if((nx != inF.getIndex1Size()+1) || (ny != inF.getIndex2Size()+1))
    {
    initialize(inF.getIndex1Size()+1,inF.getIndex2Size()+1);
    }

    for(i=0; i < nSampleX; i++)
    {
    for(j=0; j < nSampleY; j++)
    {
        in[j+ i*nSampleY] = inF(i,j);
    }}

    fftw_execute(plan);

    // Capture transform values and scale appropriately

    double scalingfactor = sqrt(LX*LY)/(2.0*((double)(nx)*(double)(ny)));
    
    for(i=0; i < nSampleX; i++)
    {
    for(j=0; j < nSampleY; j++)
    {
        outF(i,j) = out[j+ i*nSampleY]*scalingfactor;
    }}
}
    
// fftw2d_sin_inverse argument sizes:
//
// DoubleVector2d size nx-1 by ny-1
//

void fftw2d_sin_inverse(DoubleVector2d& inF, DoubleVector2d& outF)
{
    long i,j;

    if((nx != inF.getIndex1Size()+1) || (ny != inF.getIndex2Size()+1))
    {
    initialize(inF.getIndex1Size()+1,inF.getIndex2Size()+1);
    }

    for(i=0; i < nSampleX; i++)
    {
    for(j=0; j < nSampleY; j++)
    {
        in[j+ i*nSampleY] = inF(i,j);
    }}

    fftw_execute(plan);

    double scalingfactor = 1.0/(2.0*sqrt(LX*LY));

    for(i=0; i < nSampleX; i++)
    {
    for(j=0; j < nSampleY; j++)
    {
        outF(i,j) = out[j+ i*nSampleY]*scalingfactor;
    }}
}

//
// "In place sin transforms"
//
// fftw2d_sin_forward argument sizes:
//
// DoubleVector2d size nx-1 by ny-1
//
void fftw2d_sin_forward(DoubleVector2d& F)
{
    long i,j;

    if((nx != F.getIndex1Size()+1) || (ny != F.getIndex2Size()+1))
    {
    initialize(F.getIndex1Size()+1,F.getIndex2Size()+1);
    }

    for(i=0; i < nSampleX; i++)
    {
    for(j=0; j < nSampleY; j++)
    {
        in[j+ i*nSampleY] = F(i,j);
    }}

    fftw_execute(plan);

    // Capture transform values and scale appropriately

    double scalingfactor = sqrt(LX*LY)/(2.0*((double)(nx)*(double)(ny)));

    for(i=0; i < nSampleX; i++)
    {
    for(j=0; j < nSampleY; j++)
    {
        F(i,j) = out[j+ i*nSampleY]*scalingfactor;
    }}
}

// fftw2d_sin_inverse argument sizes:
//
// DoubleVector2d size nx-1 by ny-1
//

void fftw2d_sin_inverse(DoubleVector2d&  F)
{
    long i,j;

    if((nx != F.getIndex1Size()+1) || (ny != F.getIndex2Size()+1))
    {
    initialize(F.getIndex1Size()+1,F.getIndex2Size()+1);
    }

    for(i=0; i < nSampleX; i++)
    {
    for(j=0; j < nSampleY; j++)
    {
        in[j+ i*nSampleY] = F(i,j);
    }}

    fftw_execute(plan);

    double scalingfactor = 1.0/(2.0*sqrt(LX*LY));

    for(i=0; i < nSampleX; i++)
    {
    for(j=0; j < nSampleY; j++)
    {
        F(i,j) = out[j+ i*nSampleY]*scalingfactor;
    }}
}


private:

    long nx;     long ny;
    double LX; double LY;
    
    long nSampleX;
    long nSampleY;

    fftw_plan plan;

    double*  in;
    double* out;
};
}
#endif
