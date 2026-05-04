/*
 * Convolution1d.h
 *
 *  Created on: Jul 30, 2020
 *      Author: anderson
 */


#ifndef CONVOLUTION_1D_
#define CONVOLUTION_1D_

#include "GridFunctionNd/SCC_GridFunction1d.h"
#include "FFTW3_InterfaceNd/SCC_fftw3_1d.h"
#include "FFTW3_InterfaceNd/SCC_FFT_Nvalues.h"
//
// A class to compute the convolution sum of two grid functions
// F and G of the form
//
//              [k = M/2]
// (F o G)_j = Sum              F(x_j - x_k)*G(x_k)*hx
//              [k = -M/2]
//
// Here F is the "target function" defined on a domain of arbitrary size with mesh
// width hx and the kernel G is specified on a domain [-(M/2)*hx,(M/2)*hx]
// discretized with an even number of panels M.
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
class Convolution1d
{
    SCC::FFT_Nvalues fftNvalues;

    public:

    Convolution1d()
    {
    initialize();
    }

    Convolution1d(const Convolution1d& C)
    {
    if(C.panelCountX == 0) {initialize(); return;}
    initialize(C);
    }

    void initialize()
    {
    xPanelsG = 0;
    xPanelsF = 0;
    hx       = 0;
    panelCountX = 0;

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

    void initialize(const Convolution1d& C)
    {
    xPanelsG    = C.xPanelsG;
    xPanelsF    = C.xPanelsF;
    hx          = C.hx;
    panelCountX = C.panelCountX;

    Fext.initialize(C.Fext);

    Gwrap.initialize(C.Gwrap);

    DFT.initialize(panelCountX);

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
    void initialize(SCC::GridFunction1d& F, SCC::GridFunction1d& G)
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

    bool isEqualStructure(SCC::GridFunction1d& F, SCC::GridFunction1d& G)
    {
    long sizeDiff = 0;

    sizeDiff += std::abs(xPanelsF - F.getXpanelCount());
    sizeDiff += std::abs(xPanelsG - G.getXpanelCount());

    if(std::abs(hx - F.getHx()) > 1.0e-12) sizeDiff += 1;
    if(std::abs(hx - G.getHx()) > 1.0e-12) sizeDiff += 1;

    if(sizeDiff > 0) return false;
    return true;
    }

    void allocateData(SCC::GridFunction1d& F, SCC::GridFunction1d& G)
    {

    // Quick return if data sizes don't need to be changed

    if(isEqualStructure(F,G)) return;

    xPanelsG = G.getXpanelCount();
    xPanelsF = F.getXpanelCount();
    hx       = F.getHx();

    if(xPanelsG%2 != 0)
    {
    throw std::runtime_error("Convolution1d Error:\nConvolution kernel must be defined on an even number of panels.\n");
    }
    if(std::abs(hx - G.getHx()) > 1.0e-12)
    {
    throw std::runtime_error("Convolution1d Error:\nMesh widths (hx) of arguments not equal.\n");
    }

    fftNvalues.setMaximalPrimeFactor(5);
    panelCountX = xPanelsF + xPanelsG/2 + 1;
    panelCountX = fftNvalues.getFFT_N(panelCountX);

    double xStarMin =  F.getXmin();
    double xStarMax =  xStarMin + panelCountX*hx;

    Fext.initialize(panelCountX,xStarMin,xStarMax);
    Gwrap.initialize(panelCountX,xStarMin,xStarMax);

    DFT.initialize(panelCountX);

    dataReal.initialize(panelCountX,xStarMin,xStarMax);
    dataImag.initialize(panelCountX,xStarMin,xStarMax);

	transDataReal_F.initialize(panelCountX);
	transDataImag_F.initialize(panelCountX);

	transDataReal_G.initialize(panelCountX);
	transDataImag_G.initialize(panelCountX);

	transDataReal.initialize(panelCountX);
	transDataImag.initialize(panelCountX);
    }

	void FconvG(SCC::GridFunction1d& F, SCC::GridFunction1d& G)
	{

	allocateData(F,G);

    // Normalize the domain for the transforms

    double xStarMin =  F.getXmin();
    double xStarMax =  xStarMin + panelCountX*hx;

    // Embed the function to be convolved in the larger domain (Fext)

    Fext.setToValue(0.0);
    for(long i = 0; i <= xPanelsF; i++)
    {
    	Fext(i) = F(i);
    }

    //
    // Store G in wrap-around order
    //

    Gwrap.setToValue(0.0);

    long iInd;

    for(long i = -xPanelsG/2; i <= xPanelsG/2; i++)
    {
    iInd = (i < 0) ? i + panelCountX : i;
    Gwrap(iInd) = G(i + xPanelsG/2);
    }

	dataImag.setToValue(0.0);

    // Compute forward transform of the target function

	dataReal = Fext;
    DFT.fftw1d_forward(dataReal,dataImag,transDataReal_F,transDataImag_F);

    // Compute the forward transform of the mollifier

	dataReal = Gwrap;
    DFT.fftw1d_forward(dataReal,dataImag,transDataReal_G,transDataImag_G);

    // multiply the transforms being mindful of complex multiplication

    transDataReal = transDataReal_G*transDataReal_F - transDataImag_G*transDataImag_F;
    transDataImag = transDataReal_G*transDataImag_F + transDataImag_G*transDataReal_F;

    // Scaling factor because of internal scaling of fftw interfaces

    transDataReal *= sqrt(xStarMax-xStarMin);
    transDataImag *= sqrt(xStarMax-xStarMin);

    DFT.fftw1d_inverse(transDataReal,transDataImag,dataReal,dataImag);

    // Capture result

    for(long i = 0; i <= xPanelsF; i++)
    {
    	F(i) = dataReal(i);
    }

    }

    long xPanelsG;
    long xPanelsF;
    double hx;
    long panelCountX;

    SCC::fftw3_1d DFT;

    SCC::GridFunction1d dataReal;
    SCC::GridFunction1d dataImag;

	SCC::DoubleVector1d	transDataReal;
	SCC::DoubleVector1d	transDataImag;

	SCC::DoubleVector1d	transDataReal_F;
	SCC::DoubleVector1d	transDataImag_F;

	SCC::DoubleVector1d	transDataReal_G;
	SCC::DoubleVector1d	transDataImag_G;

    SCC::GridFunction1d Fext;
    SCC::GridFunction1d Gwrap;
};
}
#endif
