#include "fftw3.h"

#include "../DoubleVectorNd/SCC_DoubleVector3d.h"    // Array class header
#include "../GridFunctionNd/SCC_GridFunction3d.h"    // Grid function class header

#include <exception>


#ifndef SCC_FFTW3_SIN_3D_
#define SCC_FFTW3_SIN_3D_

//
//####################################################################
//          SCC_fftw3_sin3d.h : Version Nov. 12, 2015
//####################################################################
//
/**<p>
    This class provides forward and inverse sin transforms of three dimensional
    data stored in GridFunction3d instances or data stored in
    DoubleVector3d instances. Normalization and scaling is performed so that the
    transforms computed are close analogs to the mathematical sin transforms of
    functions defined over intervals of size LX x LY x LZ.

    The transforms are implemented using the r2r_3d(...) routines from the
    FFTW3 Fast Fourier Transform library (http://www.fftw.org/). The
    FFTW3 headers and library files must be available to create any executable
    that uses this class.

    This transform interface is designed be used with the data values of
    discrete functions with homogeneous boundary values associated with
    a uniform grid discretization consisting of nx grid panels by ny grid panels
    by nx grid panels spanning a domain of size LX x LY x LZ.

    When DoubleVector3d's are used for input and/or output their size
    is (nx-1) x (ny-1) x (nz-1); the number of interior data values associated with
    such discretizations.

    For consistency with other FFTW_InterfaceNd classes, in the
    construction or initialization of class instances, the number
    of panels nx, ny and nz of a discretization is specified.

    For computational efficiency nx, ny and nz should be chosen to be a
    product of small primes less than or equal to 13  (if possible).
    
    ===============================================================

    What's computed assuming the input data is specified as a
    GridFunction3d instance initialized with the parameters
    (nx, xMin, xMax, ny, yMin, yMax, nz, zMin, zMax):

    If (kx,ky,kz) is a wave number vector with kx, ky, kz in the ranges

    1 <= kx <=  (nx-1) , 1 <= ky <= (ny-1) and 1 <= kz <= (nz-1)
    
    then a GridFunction3d instance whose values consist of

    d  = sqrt(2/LX)*sin(kx*pi*(x_i-xMin)/LX)
        *sqrt(2/LY)*sin(ky*pi*(y_j-yMin)/LY)
        *sqrt(2/LZ)*sin(ky*pi*(z_k-zMin)/LZ)

    evaluated the points

    x_i = xMin + i*hx, i = 0,1,...,nx-1, hx = LX/nx
    y_j = yMin + j*hy, j = 0,1,...,ny-1, hy = LY/ny
    z_k = zMin + k*hz, k = 0,1,...,nz-1, hz = LZ/nz

    gets mapped by the forward transform to a DoubleVector3d of size
    (nx-1)x*(ny-1)x(nz-1), f_hat, where

    f_hat(kx-1, ky-1, kz-1) = 1

    and all other values are zero.

    Similarly the inverse transform maps the DoubleVector3d of size
    (nx-1)x*(ny-1)x(nz-1), f_hat, that is zero value except for a value of 1 in
    (kx-1,ky-1,kz-1) place to the GridFunction3d instance with values
    corresponding to

    d  = sqrt(2/LX)*sin(kx*pi*(x_i-xMin)/LX)
        *sqrt(2/LY)*sin(ky*pi*(y_j-yMin)/LY)
        *sqrt(2/LZ)*sin(ky*pi*(z_k-zMin)/LZ)

    evaluated the points

    x_i = xMin + i*hx, i = 0,1,...,nx-1, hx = LX/nx
    y_j = yMin + j*hy, j = 0,1,...,ny-1, hy = LY/ny
    z_k = zMin + k*hz, k = 0,1,...,nz-1, hz = LZ/nz

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
    for(long kz = 1; kz <= nz-1; kz++)
    {
        sCoeff  = f_hat(kx-1,ky-1,kz-1);  // (kx,ky,kz)'th  sin coefficient is
                *                         // the (kx-1,ky-1,kz-1) entry of the transform
                *
                *
    }}}
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
// Chris Anderson (C) UCLA                             Nov. 12, 2015
//####################################################################
//
namespace SCC
{
class fftw3_sin3d
{
public:

fftw3_sin3d()
{
    plan = nullptr;
    in   = nullptr;
    out  = nullptr;

    nx  = 0;
    ny  = 0;
    nz  = 0;

    LX = 1.0;
    LY = 1.0;
    LZ = 1.0;

    nSampleX = 0;
    nSampleY = 0;
    nSampleZ = 0;
}

fftw3_sin3d(long nx, long ny, long nz, double LX = 1.0, double LY = 1.0, double LZ = 1.0)
{
    this->nx    = 0;
    this->ny    = 0;
    this->nz    = 0;

    in          = nullptr;
    out         = nullptr;
    plan        = nullptr;

    initialize(nx,ny,nz,LX,LY,LZ);
}

virtual ~fftw3_sin3d()
{
    if(plan != nullptr)
    {fftw_destroy_plan(plan); /*fftw_cleanup();*/}

    if(in  != nullptr) fftw_free(in);
    if(out != nullptr) fftw_free(out);
}

void initialize()
{
    if(plan != nullptr)
    {fftw_destroy_plan(plan);}
    
    if(in   != nullptr) fftw_free(in);
    if(out  != nullptr) fftw_free(out);

    plan = nullptr;
    in   = nullptr;
    out  = nullptr;

    nx  = 0;
    ny  = 0;
    nz  = 0;

    LX = 1.0;
    LY = 1.0;
    LZ = 1.0;

    nSampleX = 0;
    nSampleY = 0;
    nSampleZ = 0;
}

void initialize(long nx, long ny, long nz, double LX = 1.0, double LY = 1.0, double LZ = 1.0)
{
    if((this->nx != nx)||(this->ny != ny)||(this->nz != nz))
    {
    this->nx   = nx;
    this->ny   = ny;
    this->nz   = nz;
    nSampleX   = nx-1;
    nSampleY   = ny-1;
    nSampleZ   = nz-1;

    if(plan != 0)
    {fftw_destroy_plan(plan);}

    if(in  != 0) fftw_free(in);
    if(out != 0) fftw_free(out);

    in  = (double*) fftw_malloc(sizeof(double) * nSampleX*nSampleY*nSampleZ);
    out = (double*) fftw_malloc(sizeof(double) * nSampleX*nSampleY*nSampleZ);

    plan = fftw_plan_r2r_3d(nSampleX, nSampleY, nSampleZ, in, out,FFTW_RODFT00,
           FFTW_RODFT00, FFTW_RODFT00,FFTW_ESTIMATE);

    if(plan == nullptr)
    {
    throw std::runtime_error("\nXXX Error : required fftw_r2r_3d function not available \nXXX in FFTW library used (likely MKL) ");
    }

    // Efficiency: From the documentation for RODFT00 transform, nSampleX+1, nSampleY+1, nSampleZ+1
    // should be a product of small primes ==> nx, ny and nz should be product of small primes.
    }

    this->LX = LX;
    this->LY = LY;
    this->LZ = LZ;
}


// replan() is used to create new plans for multi-threaded instances when
// the number of threads allocated to FFTW3 is changed.

void replan()
{
   if(plan != nullptr) {fftw_destroy_plan(plan);}

   if((in  == nullptr) || (out == nullptr))
   {
   throw std::runtime_error("\nXXX Error : calling replan() before SCC::fft3_sin3d instance initialized.\n ");
   }

   plan = fftw_plan_r2r_3d(nSampleX, nSampleY, nSampleZ, in, out,FFTW_RODFT00, FFTW_RODFT00, FFTW_RODFT00,FFTW_ESTIMATE);

   if(plan == nullptr)
   {
    throw std::runtime_error("\nXXX Error : required fftw_r2r_3d function not available \nXXX in FFTW library used (likely MKL) ");
   }
}




// fftw3d_sin_forward argument sizes:
//
// DoubleVector3d size (nx-1) x  (ny-1) x (nz-1)
//
// This operator ignores the perimeter values of the input GridFunction3d

void fftw3d_sin_forward(GridFunction3d& inF, DoubleVector3d& outF)
{
    this->LX = inF.getXmax() - inF.getXmin();
     this->LY = inF.getYmax() - inF.getYmin();
     this->LZ = inF.getZmax() - inF.getZmin();

    long i,j,k;

    if((nx != inF.getXpanelCount()) || (ny != inF.getYpanelCount())|| (nz != inF.getZpanelCount()))
    {
    initialize(inF.getXpanelCount(),inF.getYpanelCount(),inF.getZpanelCount(),LX,LY,LZ);
    }

    for(i=0; i < nSampleX; i++)
    {
    for(j=0; j < nSampleY; j++)
    {
    for(k=0; k < nSampleZ; k++)
    {
        in[k + nSampleZ*(j+ i*nSampleY)] = inF(i+1,j+1,k+1);
    }}}

    fftw_execute(plan);

    // Capture transform values and scale appropriately

    double scalingfactor =  sqrt(LX*LY*LZ)/(2.0*sqrt(2.0)*((double)(nx)*(double)(ny)*(double)(nz)));

    for(i=0; i < nSampleX; i++)
    {
    for(j=0; j < nSampleY; j++)
    {
    for(k=0; k < nSampleZ; k++)
    {
        outF(i,j,k) = out[k + nSampleZ*(j+ i*nSampleY)]*scalingfactor;
    }}}
}

// fftw3d_sin_inverse argument sizes:
//
// DoubleVector3d size (nx-1) x  (ny-1) x (nz-1)
//

void fftw3d_sin_inverse(DoubleVector3d& inF, GridFunction3d& outF)
{
    long i,j,k;

    this->LX = outF.getXmax() - outF.getXmin();
    this->LY = outF.getYmax() - outF.getYmin();
    this->LZ = outF.getZmax() - outF.getZmin();

    if((nx != inF.getIndex1Size()+1) || (ny != inF.getIndex2Size()+1)|| (nz != inF.getIndex3Size()+1))
    {
    initialize(inF.getIndex1Size()+1,inF.getIndex2Size()+1,inF.getIndex3Size()+1,LX,LY,LZ);
    }

    for(i=0; i < nSampleX; i++)
    {
    for(j=0; j < nSampleY; j++)
    {
    for(k=0; k < nSampleZ; k++)
    {
        in[k + nSampleZ*(j+ i*nSampleY)] = inF(i,j,k);
    }}}

    fftw_execute(plan);

    double scalingfactor = 1.0/(2.0*sqrt(2.0*LX*LY*LZ));

    for(i=0; i < nSampleX; i++)
    {
    for(j=0; j < nSampleY; j++)
    {
    for(k=0; k < nSampleZ; k++)
    {
        outF(i+1,j+1,k+1) = out[k + nSampleZ*(j+ i*nSampleY)]*scalingfactor;
    }}}

    outF.setBoundaryValues(0.0);
}


    
// fftw3d_sin_forward argument sizes:
//
// DoubleVector3d size (nx-1) x  (ny-1) x (nz-1)
//

void fftw3d_sin_forward(DoubleVector3d& F)
{
    long i,j,k;

    if((nx != F.getIndex1Size()+1) || (ny != F.getIndex2Size()+1)|| (nz != F.getIndex3Size()+1))
    {
    initialize(F.getIndex1Size()+1,F.getIndex2Size()-1,F.getIndex3Size()+1);
    }

    for(i=0; i < nSampleX; i++)
    {
    for(j=0; j < nSampleY; j++)
    {
    for(k=0; k < nSampleZ; k++)
    {
        in[k + nSampleZ*(j+ i*nSampleY)] = F(i,j,k);
    }}}

    fftw_execute(plan);

    // Capture transform values and scale appropriately

    double scalingfactor = sqrt(LX*LY*LZ)/(2.0*sqrt(2.0)*((double)(nx)*(double)(ny)*(double)(nz)));

    for(i=0; i < nSampleX; i++)
    {
    for(j=0; j < nSampleY; j++)
    {
    for(k=0; k < nSampleZ; k++)
    {
        F(i,j,k) = out[k + nSampleZ*(j+ i*nSampleY)]*scalingfactor;
    }}}
}

// fftw3d_sin_forward argument sizes:
//
// DoubleVector3d size (nx-1) x  (ny-1) x (nz-1)
//

void fftw3d_sin_forward(DoubleVector3d& inF, DoubleVector3d& outF)
{
    long i,j,k;

    if((nx != inF.getIndex1Size()+1) || (ny != inF.getIndex2Size()+1)|| (nz != inF.getIndex3Size()+1))
    {
    initialize(inF.getIndex1Size()+1,inF.getIndex2Size()+1,inF.getIndex3Size()+1);
    }

    for(i=0; i < nSampleX; i++)
    {
    for(j=0; j < nSampleY; j++)
    {
    for(k=0; k < nSampleZ; k++)
    {
        in[k + nSampleZ*(j+ i*nSampleY)] = inF(i,j,k);
    }}}

    fftw_execute(plan);

    // Capture transform values and scale appropriately
    
    double scalingfactor =  sqrt(LX*LY*LZ)/(2.0*sqrt(2.0)*((double)(nx)*(double)(ny)*(double)(nz)));

    for(i=0; i < nSampleX; i++)
    {
    for(j=0; j < nSampleY; j++)
    {
    for(k=0; k < nSampleZ; k++)
    {
        outF(i,j,k) = out[k + nSampleZ*(j+ i*nSampleY)]*scalingfactor;
    }}}
}
    
// fftw3d_sin_inverse argument sizes:
//
// DoubleVector3d size (nx-1) x  (ny-1) x (nz-1)
//

void fftw3d_sin_inverse(DoubleVector3d& inF, DoubleVector3d& outF)
{
    long i,j,k;

    if((nx != inF.getIndex1Size()+1) || (ny != inF.getIndex2Size()+1)|| (nz != inF.getIndex3Size()+1))
    {
    initialize(inF.getIndex1Size()+1,inF.getIndex2Size()+1,inF.getIndex3Size()+1);
    }

    for(i=0; i < nSampleX; i++)
    {
    for(j=0; j < nSampleY; j++)
    {
    for(k=0; k < nSampleZ; k++)
    {
        in[k + nSampleZ*(j+ i*nSampleY)] = inF(i,j,k);
    }}}

    fftw_execute(plan);

    double scalingfactor = 1.0/(2.0*sqrt(2.0*LX*LY*LZ));

    for(i=0; i < nSampleX; i++)
    {
    for(j=0; j < nSampleY; j++)
    {
    for(k=0; k < nSampleZ; k++)
    {
        outF(i,j,k) = out[k + nSampleZ*(j+ i*nSampleY)]*scalingfactor;
    }}}
}

// fftw3d_sin_inverse argument sizes:
//
// DoubleVector3d size (nx-1) x  (ny-1) x (nz-1)
//
    
void fftw3d_sin_inverse(DoubleVector3d&  F)
{
    long i,j,k;

    if((nx != F.getIndex1Size()+1) || (ny != F.getIndex2Size()+1)|| (nz != F.getIndex3Size()+1))
    {
    initialize(F.getIndex1Size()+1,F.getIndex2Size()+1,F.getIndex3Size()+1);
    }

    for(i=0; i < nSampleX; i++)
    {
    for(j=0; j < nSampleY; j++)
    {
    for(k=0; k < nSampleZ; k++)
    {
        in[k + nSampleZ*(j+ i*nSampleY)] = F(i,j,k);
    }}}

    fftw_execute(plan);

    double scalingfactor = 1.0/(2.0*sqrt(2.0*LX*LY*LZ));

    for(i=0; i < nSampleX; i++)
    {
    for(j=0; j < nSampleY; j++)
    {
    for(k=0; k < nSampleZ; k++)
    {
        F(i,j,k) = out[k + nSampleZ*(j+ i*nSampleY)]*scalingfactor;
    }}}
}


private:


    long nx;   long ny;     long nz;
    double LX; double LY; double LZ;
    
    long nSampleX;
    long nSampleY;
    long nSampleZ;

    fftw_plan plan;

    double*  in;
    double* out;
};
}
#endif
