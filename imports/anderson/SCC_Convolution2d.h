/*
 * Convolution2d.h
 *
 *  Created on: Jul 31, 2020
 *      Author: anderson
 */

#ifndef CONVOLUTION_2D_
#define CONVOLUTION_2D_

#include "GridFunctionNd/SCC_GridFunction2d.h"
#include "FFTW3_InterfaceNd/SCC_fftw3_2d.h"
#include "FFTW3_InterfaceNd/SCC_FFT_Nvalues.h"
//
// A class to compute the convolution sum of two grid functions
// F and G of the form
//
//                   [q = MY/2]
//                   [p = MX/2]
// (F o G)_(i,j)  = Sum              F(x_i - x_p, y_j  - y_q)*G(x_p,y_q)*hx*hy
//                   [p = -MX/2]
//                   [q = -MY/2]
//
// Here F is the "target function" defined on a domain of arbitrary size with mesh
// width hx and the kernel G is specified on a domain [-(MX/2)*hx,(MX/2)*hx] X [-(MY/2)*hy,(MY/2)*hy]
// discretized with an even number of panels (MX,MY) in each direction. .
//
/*
#############################################################################
#
# Copyright 2020 Chris Anderson
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
class Convolution2d
{
    SCC::FFT_Nvalues fftNvalues;

    public:

    Convolution2d()
    {
    initialize();
    }

    Convolution2d(const Convolution2d& C)
    {
    if((C.panelCountX == 0)||(C.panelCountY == 0)) {initialize(); return;}
    initialize(C);
    }

    void initialize()
    {
    xPanelsG = 0;
    xPanelsF = 0;
    hx       = 0;
    panelCountX = 0;

    yPanelsG = 0;
    yPanelsF = 0;
    hy       = 0;
    panelCountY = 0;
    Fext.initialize();

    Gwrap.initialize();

    dataReal.initialize();
    dataImag.initialize();

	transDataReal_F.initialize();
	transDataImag_F.initialize();

	transDataReal_G.initialize();
	transDataImag_G.initialize();

	transDataReal.initialize();
	transDataImag.initialize();
    }

    void initialize(const Convolution2d& C)
    {
    xPanelsG    = C.xPanelsG;
    xPanelsF    = C.xPanelsF;
    hx          = C.hx;
    panelCountX = C.panelCountX;

    yPanelsG    = C.yPanelsG;
    yPanelsF    = C.yPanelsF;
    hy          = C.hy;
    panelCountY = C.panelCountY;

    Fext.initialize(C.Fext);

    Gwrap.initialize(C.Gwrap);

    DFT.initialize(panelCountX,panelCountY);

    dataReal.initialize(C.dataReal);
    dataImag.initialize(C.dataImag);

	transDataReal_F.initialize(C.transDataReal_F);
	transDataImag_F.initialize(C.transDataImag_F);

	transDataReal_G.initialize(C.transDataReal_G);
	transDataImag_G.initialize(C.transDataImag_G);

	transDataReal.initialize(C.transDataReal);
	transDataImag.initialize(C.transDataImag);
    }

    // Initializes data structures based upon structure of target F
    // and kernal G
    ///
    void initialize(SCC::GridFunction2d& F, SCC::GridFunction2d& G)
    {
    	allocateData(F,G);
    }

    void destroyFFTWplans()
    {
    DFT.initialize();
    }

// replan() is used to create new plans for multi-threaded instances when
// the number of threads allocated to FFTW3 is changed.

    void replan()
    {
    DFT.replan();
    }


    // Returns true if the input arguments have data structure parameters
    // identical to the cached version (so no need to reallocate data).

    bool isEqualStructure(SCC::GridFunction2d& F, SCC::GridFunction2d& G)
    {
    long sizeDiff = 0;

    sizeDiff += std::abs(xPanelsF - F.getXpanelCount());
    sizeDiff += std::abs(yPanelsF - F.getYpanelCount());
    sizeDiff += std::abs(xPanelsG - G.getXpanelCount());
    sizeDiff += std::abs(yPanelsG - G.getYpanelCount());

    if(std::abs(hx - F.getHx()) > 1.0e-12) sizeDiff += 1;
    if(std::abs(hy - F.getHy()) > 1.0e-12) sizeDiff += 1;
    if(std::abs(hx - G.getHx()) > 1.0e-12) sizeDiff += 1;
    if(std::abs(hy - G.getHy()) > 1.0e-12) sizeDiff += 1;

    if(sizeDiff > 0) return false;
    return true;
    }


    void allocateData(SCC::GridFunction2d& F, SCC::GridFunction2d& G)
    {

    // Quick return if data sizes don't need to be changed

    if(isEqualStructure(F,G)) return;

    xPanelsG = G.getXpanelCount();
    xPanelsF = F.getXpanelCount();
    hx       = F.getHx();

    yPanelsG = G.getYpanelCount();
    yPanelsF = F.getYpanelCount();
    hy       = F.getHy();

    if((xPanelsG%2 != 0)||(yPanelsG%2 != 0))
    {
    throw std::runtime_error("Convolution1d Error:\nConvolution kernel must be defined on an even number of panels.\n");
    }
    if((std::abs(hx - G.getHx()) > 1.0e-12)||(std::abs(hy - G.getHy()) > 1.0e-12 ))
    {
    throw std::runtime_error("Convolution2d Error:\nMesh widths (hx,hy) of arguments not equal.\n");
    }

    fftNvalues.setMaximalPrimeFactor(7);
    
    panelCountX = xPanelsF + xPanelsG/2 + 1;
    panelCountY = yPanelsF + yPanelsG/2 + 1;
    panelCountX = fftNvalues.getFFT_N(panelCountX);
    panelCountY = fftNvalues.getFFT_N(panelCountY);

    double xStarMin =  F.getXmin();
    double xStarMax =  xStarMin + panelCountX*hx;

    double yStarMin =  F.getYmin();
    double yStarMax =  yStarMin + panelCountY*hy;

    Fext.initialize(panelCountX,xStarMin,xStarMax,panelCountY,yStarMin,yStarMax);
    Gwrap.initialize(panelCountX,xStarMin,xStarMax,panelCountY,yStarMin,yStarMax);

    DFT.initialize(panelCountX,panelCountY);

    dataReal.initialize(panelCountX,xStarMin,xStarMax,panelCountY,yStarMin,yStarMax);
    dataImag.initialize(panelCountX,xStarMin,xStarMax,panelCountY,yStarMin,yStarMax);

	transDataReal_F.initialize(panelCountX,panelCountY);
	transDataImag_F.initialize(panelCountX,panelCountY);

	transDataReal_G.initialize(panelCountX,panelCountY);
	transDataImag_G.initialize(panelCountX,panelCountY);

	transDataReal.initialize(panelCountX,panelCountY);
	transDataImag.initialize(panelCountX,panelCountY);
    }

	void FconvG(SCC::GridFunction2d& F, SCC::GridFunction2d& G)
	{

	allocateData(F,G);

    // Normalize the domain for the transforms

    double xStarMin =  F.getXmin();
    double xStarMax =  xStarMin + panelCountX*hx;

    double yStarMin =  F.getYmin();
    double yStarMax =  yStarMin + panelCountY*hy;

    // Embed the function to be convolved in the larger domain (Fext)

    Fext.setToValue(0.0);
    for(long i = 0; i <= xPanelsF; i++)
    {
    for(long j = 0; j <= yPanelsF; j++)
    {
    	Fext(i,j) = F(i,j);
    }}

    //
    // Store G in wrap-around order
    //

    Gwrap.setToValue(0.0);

    long iInd; long jInd;

    for(long i = -xPanelsG/2; i <= xPanelsG/2; i++)
    {
    iInd = (i < 0) ? i + panelCountX : i;
    for(long j = -yPanelsG/2; j <= yPanelsG/2; j++)
    {
    jInd = (j < 0) ? j + panelCountY : j;
    Gwrap(iInd,jInd) = G(i + xPanelsG/2,j+yPanelsG/2);
    }}

	dataImag.setToValue(0.0);

    // Compute forward transform of the target function

	dataReal = Fext;
    DFT.fftw2d_forward(dataReal,dataImag,transDataReal_F,transDataImag_F);

    // Compute the forward transform of the mollifier

	dataReal = Gwrap;
    DFT.fftw2d_forward(dataReal,dataImag,transDataReal_G,transDataImag_G);

    // multiply the transforms being mindful of complex multiplication

    transDataReal = transDataReal_G*transDataReal_F - transDataImag_G*transDataImag_F;
    transDataImag = transDataReal_G*transDataImag_F + transDataImag_G*transDataReal_F;

    // Scaling factor because of internal scaling of fftw interfaces

    transDataReal *= sqrt(xStarMax-xStarMin)*sqrt(yStarMax-yStarMin);
    transDataImag *= sqrt(xStarMax-xStarMin)*sqrt(yStarMax-yStarMin);

    DFT.fftw2d_inverse(transDataReal,transDataImag,dataReal,dataImag);

    // Capture result

    for(long i = 0; i <= xPanelsF; i++)
    {
    for(long j = 0; j <= yPanelsF; j++)
    {
    	F(i,j) = dataReal(i,j);
    }}

    }

    long xPanelsG;
    long xPanelsF;
    double hx;
    long panelCountX;

    long yPanelsG;
    long yPanelsF;
    double hy;
    long panelCountY;

    SCC::fftw3_2d DFT;

    SCC::GridFunction2d dataReal;
    SCC::GridFunction2d dataImag;

	SCC::DoubleVector2d	transDataReal;
	SCC::DoubleVector2d	transDataImag;

	SCC::DoubleVector2d	transDataReal_F;
	SCC::DoubleVector2d	transDataImag_F;

	SCC::DoubleVector2d	transDataReal_G;
	SCC::DoubleVector2d	transDataImag_G;

    SCC::GridFunction2d    Fext;
    SCC::GridFunction2d   Gwrap;
};
}
#endif




