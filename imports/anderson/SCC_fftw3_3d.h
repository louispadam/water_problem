#include "fftw3.h"
#include "../DoubleVectorNd/SCC_DoubleVector3d.h"
#include "../GridFunctionNd/SCC_GridFunction3d.h"

#ifndef SCC_FFTW3_3D_
#define SCC_FFTW3_3D_
//
//####################################################################
//          fftw3_3d.h : Version  Nov. 12, 2015
//####################################################################
//
/**<p>
    This class provides forward and inverse transforms of three dimensional
    data stored in GridFunction3d instances or data stored in
    DoubleVector3d instances. Normalization and scaling is performed so that the
    transforms computed are close analogs to the mathematical transforms of
    periodic functions defined over intervals of size LX x LY x LZ.

    The transforms are implemented using the dft_3d(...) routines from the
    FFTW3 Fast Fourier Transform library (http://www.fftw.org/). The
    FFTW3 headers and library files must be available to create any executable
    that uses this class.

    This transform interface is designed be used with the data values of
    discrete  periodic functions associated with a uniform grid discretization
    consisting of nx grid panels by ny grid panels by nz grid panels
    spanning a rectangular domain of size LX x LY x LZ.

    When DoubleVector3d's are used for input and/or output their
    size is nx x ny x nz; the number of independent data values associated
    with the representation of periodic functions of such discretizations.

    For computational efficiency nx, ny and nz should be chosen to be a
    product of small primes less than or equal to 13  (if possible).

    ===============================================================

    What's computed assuming the input data is specified as a
    GridFunction3d instance initialized with the parameters
    (nx, xMin, xMax, ny, yMin, yMax, nz, zMin, zMax):

    If (kx,ky,kz) is a wave number vector with kx,ky and kz in the
    ranges

    -(nx/2) <= kx <=  (nx-1)/2
    -(ny/2) <= ky <=  (ny-1)/2,
    -(nz/2) <= kz <=  (nz-1)/2

    (the range calculated with integer division), then two GridFunction3d
    instances whose values are the real and imaginary parts of

    d  = sqrt(1/LX)*exp(2*pi*I*kx*(x-xMin)/LX)
        *sqrt(1/LY)*exp(2*pi*I*ky*(y-yMin)/LY)
        *sqrt(1/LZ)*exp(2*pi*I*kz*(z-zMin)/LZ)

    evaluated the points

    x = xMin + i*hx, i = 0,1,...,nx-1, hx = LX/nx
    y = yMin + j*hy, j = 0,1,...,ny-1, hy = LY/ny
    z = zMin + k*hz, j = 0,1,...,nz-1, hz = LZ/nz

    get mapped by the forward transform to the two DoubleVector3d instances
    of size nx x ny x nz, fhat_real and fhat_imag, where 

    fhat_real(kx + (nx/2), ky + (ny/2),kz + (nz/2)) = 1

    and all other values of fhat_real and fhat_imag are zero.

    Similarly the inverse transform maps DoubleVector3d instances fhat_real and fhat_imag
    that are zero except for a value 1 in the position (kx + (nx/2), ky + (ny/2),kz + (nz/2))
    to two GridFunction3d instances whose values consist of the real and imaginary components
    of

    d  = sqrt(1/LX)*exp(2*pi*I*kx*(x-xMin)/LX)
        *sqrt(1/LY)*exp(2*pi*I*ky*(y-yMin)/LY)
        *sqrt(1/LZ)*exp(2*pi*I*kz*(z-zMin)/LZ)

    evaluated at the points

    x = xMin + i*hx, i = 0,1,...,nx-1, hx = LX/nx
    y = yMin + j*hy, j = 0,1,...,ny-1, hy = LY/ny
    z = zMin + k*hz, j = 0,1,...,nz-1, hz = LZ/nz


    After computing the transform of the data, one often wants to
    work with the discrete Fourier coefficients. Typically this is done by
    looping over the wave numbers, and then obtaining the coefficient by
    computing the correct offset. For example, one uses code constructs
    as indicated by the fragment below.

    //
    // Loop over wave numbers
    //

    long kxIndex; long kyIndex; long kzIndex;

    for(kx = -(nx/2); kx <= (nx-1)/2; kx++)
    {
    for(ky = -(ny/2); ky <= (ny-1)/2; ky++)
    {
    for(kz = -(nz/2); kz <= (nz-1)/2; kz++)
    {
    kxIndex = kx + (nx/2);             // (kxIndex, kyIndex,kzIndex) = index of the transform
    kyIndex = ky + (ny/2);             //  coefficient for the (kx,ky,kz) mode.
    kzIndex = kz + (nz/2);             //

    rk = realCoeff(kxIndex,kyIndex,kzIndex);   //  real coefficient of the (kx,ky,kz) mode
    ik = imagCoeff(kxIndex,kyIndex,kzIndex);   //  imag coefficient of the (kx,ky,kz) mode

             *
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

// See comments at the bottom of this file for information about FFTW3
// plans when using OpenMP

namespace SCC
{
class fftw3_3d
{
public:

fftw3_3d()
{
    forwardplan = nullptr;
    inverseplan = nullptr;

    in  = nullptr;
    out = nullptr;

    nx  = 0;
    ny  = 0;
    nz  = 0;

    LX = 1.0;
    LY = 1.0;
    LZ = 1.0;
}

fftw3_3d(const fftw3_3d& DFT)
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

     initialize(DFT.nx,DFT.ny,DFT.nz,DFT.LX,DFT.LY,DFT.LZ);
}

fftw3_3d(long nx, long ny, long nz, double LX = 1.0, double LY = 1.0, double LZ = 1.0)
{
    this->nx    = 0;
    this->ny    = 0;
    this->nz    = 0;

    in          = nullptr;
    out         = nullptr;
    forwardplan = nullptr;
    inverseplan = nullptr;

    initialize(nx,ny,nz,LZ,LY,LZ);
}

virtual ~fftw3_3d()
{
    if(forwardplan != nullptr)
    { fftw_destroy_plan(forwardplan);}
    
    if(inverseplan != nullptr)
    {fftw_destroy_plan(inverseplan);}

    
    if(in  != nullptr) fftw_free(in);
    if(out != nullptr) fftw_free(out);
}

void initialize()
{
    if(forwardplan != nullptr)
    { fftw_destroy_plan(forwardplan);}
    
    if(inverseplan != nullptr)
    {fftw_destroy_plan(inverseplan);}
    
    if(in  != nullptr) fftw_free(in);
    if(out != nullptr) fftw_free(out);

    forwardplan = nullptr;
    inverseplan = nullptr;

    in  = nullptr;
    out = nullptr;

    nx  = 0;
    ny  = 0;   
    nz  = 0;

    LX = 1.0;
    LY = 1.0;
    LZ = 1.0;
}

void initialize(long nx, long ny, long nz, double LX = 1.0, double LY = 1.0, double LZ = 1.0)
{
    if((this->nx != nx)||(this->ny != ny)||(this->nz != nz))
    {
    this->nx = nx;
    this->ny = ny;
    this->nz = nz;

    if(forwardplan != nullptr)
    { fftw_destroy_plan(forwardplan);}
    
    if(inverseplan != nullptr)
    {fftw_destroy_plan(inverseplan);}

    if(in  != nullptr) fftw_free(in);
    if(out != nullptr) fftw_free(out);
    
    in  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx*ny*nz);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx*ny*nz);
    

    forwardplan = fftw_plan_dft_3d(nx, ny, nz, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    inverseplan = fftw_plan_dft_3d(nx, ny, nz, in, out, FFTW_BACKWARD,FFTW_ESTIMATE);
    
    //
    // Note: In the FFTW_BACKWARD plan the "out" argument is transformed 
    //       into the "in" argument when the plan is executed. 
    //
    }

    this->LX = LX;
    this->LY = LY;
    this->LZ = LZ;
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
   throw std::runtime_error("\nXXX Error : calling replan() before SCC::fft3_3d instance initialized.\n ");
   }

    forwardplan = fftw_plan_dft_3d(nx, ny, nz, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    inverseplan = fftw_plan_dft_3d(nx, ny, nz, in, out, FFTW_BACKWARD,FFTW_ESTIMATE);
}

// fftw3d_forward argument sizes:
//
// GridFunction3d: Initialized with nx by ny by nz panels
// DoubleVector3d: Initialized of size nx by ny by nz
//

void fftw3d_forward(GridFunction3d& inReal, GridFunction3d& inImag,
DoubleVector3d& outReal, DoubleVector3d& outImag)
{
    long i,j,k,p;
    long oldindex1, oldindex2, oldindex3;

    this->LX = inReal.getXmax() - inReal.getXmin();
    this->LY = inReal.getYmax() - inReal.getYmin();
    this->LZ = inReal.getZmax() - inReal.getZmin();


    if((nx != inReal.getXpanelCount())
    || (ny != inReal.getYpanelCount())
    || (nz != inReal.getZpanelCount()))
    {
    initialize(inReal.getXpanelCount(),
               inReal.getYpanelCount(),
               inReal.getZpanelCount(),LX,LY,LZ);
    }

    // Capture and re-order input. This extraction
    // ignores the periodically defined values at
    // the outer edges of the computational domain.

    for(i=0; i < nx; i++)
    {
    for(j=0; j < ny; j++)
    {
    for(k=0; k < nz; k++)
    {
        p =  k + nz*(j + i*ny);

        in[p][0] = inReal(i,j,k);
        in[p][1] = inImag(i,j,k);
    }}}

    fftw_execute(forwardplan);

    //reorder AND scale output

    double scalingfactor = sqrt(LX*LY*LZ)/(nx*ny*nz);


    for(p=0; p < nx*ny*nz; p++)
    {

        oldindex1 = (p/nz)/ny;
        oldindex2 = (p/nz)%ny;
        oldindex3 = p%nz;

        i = (oldindex1 + (nx/2) )%nx ;//NB: uses integer division of odd numbers
        j = (oldindex2 + (ny/2) )%ny ;//NB: uses integer division of odd numbers
        k = (oldindex3 + (nz/2) )%nz ;//NB: uses integer division of odd numbers

        outReal(i,j,k) = out[p][0]*scalingfactor;
        outImag(i,j,k) = out[p][1]*scalingfactor;
    }
}


// fftw2d_inverse argument sizes:
//
// DoubleVector2d: Initialized of size nx by ny by nz
// GridFunction2d: Initialized with nx by ny by nz panels
//


void fftw3d_inverse(DoubleVector3d& inReal, DoubleVector3d& inImag,
GridFunction3d& outReal, GridFunction3d& outImag)
{
    long i,j,k;
    long p;

    this->LX = outReal.getXmax() - outReal.getXmin();
    this->LY = outReal.getYmax() - outReal.getYmin();
    this->LZ = outReal.getZmax() - outReal.getZmin();

    long newindex1,newindex2,newindex3;

    if((nx != inReal.getIndex1Size())
    || (ny != inReal.getIndex2Size())
    || (nz != inReal.getIndex3Size()))
    {
    initialize( inReal.getIndex1Size(),
                inReal.getIndex2Size(),
                inReal.getIndex3Size(),LX,LY,LZ);
    }

    for(i=0; i < nx; i++)
    {
    for(j=0; j < ny; j++)
    {
    for(k=0; k < nz; k++)
    {
        newindex1 = (i - (nx/2) + nx )%nx ;//NB: uses integer division of odd numbers
        newindex2 = (j - (ny/2) + ny )%ny ;//NB: uses integer division of odd numbers
        newindex3 = (k - (nz/2) + nz )%nz ;//NB: uses integer division of odd numbers

        p =  newindex3 + nz*(newindex2 + newindex1*ny);

        in[p][0] = inReal( i, j, k);
        in[p][1] = inImag( i, j, k);
    }}}

    //transform

    fftw_execute(inverseplan);

    double scalefactor = 1.0/sqrt(LX*LY*LZ);

    for(p=0; p < nx*ny*nz; p++)
    {
        //  p =  k + nz*(j + i*ny) = k + nz*j + i*nz*ny
        //  =>  p/nz =  j + i*ny
        //

        i = (p/nz)/ny;
        j = (p/nz)%ny;
        k = p%nz;

        outReal(i,j,k) = out[p][0]*scalefactor;
        outImag(i,j,k) = out[p][1]*scalefactor;
    }

    // Enforce periodicity of output GridFunction3d's

    outReal.enforcePeriodicity();
    outImag.enforcePeriodicity();
}
// fftw3d_forward argument sizes:
//
// DoubleVector3d: Initialized of size nx by ny by nz
//

void fftw3d_forward(DoubleVector3d& inReal, DoubleVector3d& inImag,
DoubleVector3d& outReal, DoubleVector3d& outImag)
{
    long i,j,k,p;
    long oldindex1, oldindex2, oldindex3;

    if((nx != inReal.getIndex1Size())
    || (ny != inReal.getIndex2Size())
    || (nz != inReal.getIndex3Size()))
    {
    initialize(inReal.getIndex1Size(),
               inReal.getIndex2Size(),
               inReal.getIndex3Size());
    }

    //reorder input

    for(i=0; i < nx; i++)
    {
    for(j=0; j < ny; j++)
    {
    for(k=0; k < nz; k++)
    {
        p =  k + nz*(j + i*ny);

        in[p][0] = inReal(i,j,k);
        in[p][1] = inImag(i,j,k);
    }}}

    fftw_execute(forwardplan);

    //reorder AND scale output

    double scalingfactor = sqrt(LX*LY*LZ)/(nx*ny*nz);

    for(p=0; p < nx*ny*nz; p++)
    {

        oldindex1 = (p/nz)/ny;
        oldindex2 = (p/nz)%ny;
        oldindex3 = p%nz;

        i = (oldindex1 + (nx/2) )%nx ;//NB: uses integer division of odd numbers
        j = (oldindex2 + (ny/2) )%ny ;//NB: uses integer division of odd numbers
        k = (oldindex3 + (nz/2) )%nz ;//NB: uses integer division of odd numbers

        outReal(i,j,k) = out[p][0]*scalingfactor;
        outImag(i,j,k) = out[p][1]*scalingfactor;
    }
}

// fftw3d_inverse argument sizes:
//
// DoubleVector3d: Initialized of size nx by ny by nz
//
void fftw3d_inverse(DoubleVector3d& inReal, DoubleVector3d& inImag,
DoubleVector3d& outReal, DoubleVector3d& outImag)
{
    long i,j,k;
    long p;

    long newindex1,newindex2,newindex3;

    if((nx != inReal.getIndex1Size())
    || (ny != inReal.getIndex2Size())
    || (nz != inReal.getIndex3Size()))
    {
    initialize( inReal.getIndex1Size(),
                inReal.getIndex2Size(),
                inReal.getIndex3Size());
    }

    for(i=0; i < nx; i++)
    {
    for(j=0; j < ny; j++)
    {
    for(k=0; k < nz; k++)
    {
        newindex1 = (i - (nx/2) + nx )%nx ;//NB: uses integer division of odd numbers
        newindex2 = (j - (ny/2) + ny )%ny ;//NB: uses integer division of odd numbers
        newindex3 = (k - (nz/2) + nz )%nz ;//NB: uses integer division of odd numbers

        p =  newindex3 + nz*(newindex2 + newindex1*ny);

        in[p][0] = inReal( i, j, k);
        in[p][1] = inImag( i, j, k);
    }}}

    //transform

    fftw_execute(inverseplan);

    double scaleingfactor = 1.0/sqrt(LX*LY*LZ);

    for(p=0; p < nx*ny*nz; p++)
    {
        //  p =  k + nz*(j + i*ny) = k + nz*j + i*nz*ny
        //  =>  p/nz =  j + i*ny
        //

        i = (p/nz)/ny;
        j = (p/nz)%ny;
        k = p%nz;

        outReal(i,j,k) = out[p][0]*scaleingfactor;
        outImag(i,j,k) = out[p][1]*scaleingfactor;
    }
}



    long   nx;   long ny;   long nz;
    double LX; double LY; double LZ;

    fftw_plan forwardplan;
    fftw_plan inverseplan;

    fftw_complex*  in;
    fftw_complex* out;
};
}
#endif


