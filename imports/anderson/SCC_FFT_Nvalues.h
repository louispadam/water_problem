//
// FFT_Nvalues.h 
//
// The member function getFFT_N(long nOriginal) of this class is used to determine
// an integral value of n that is larger than nOriginal and is a product of small
// primes. Currently the default maximal prime used is 13, although this maximal
// value can be reset with the member function setMaximalPrimeFactor(long factor).
//
// The value of 13 is the maximal prime occurring in a factorization of the number
// of data values whose fft's are being computed using FFTW-3 routines. If one
// is using FFTW-2, then the maximal prime should be set to 7.
//
//
// Author: Chris Anderson  
// (C) UCLA 2012 
//
// Updated:
// November 10, 2015: Added capability to set maximal prime factor (Chris Anderson)
//
/*
#############################################################################
#
# Copyright 2012-2017 Chris Anderson
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

#include <cstdio>
#include <cmath>
#include <cstdlib>
 

#ifndef SCC_FFT_VALUES_
#define SCC_FFT_VALUES_
#define DEFAULT_MAXIMAL_PRIME_INDEX_ 6
namespace SCC
{
class FFT_Nvalues
{

public:

	FFT_Nvalues()
	{
	primeIndex = DEFAULT_MAXIMAL_PRIME_INDEX_;
	}

	//
	// Returns a value of n >= nOriginal that is a product
	// of primes that are less than the prime with index primeIndex
	// (2 has primeIndex 1).
	//

	long getFFT_N(long nOriginal)
	{
	if(nOriginal == 1) return nOriginal;

	bool  checkPrime = checkLargePrimeFactors(nOriginal);
    long           n = nOriginal;

    if(checkPrime)
    {
    	while(checkPrime)
    	{
    		n++;
    		checkPrime = checkLargePrimeFactors(n);
    	}
    }
    return n;
	}
//
//  This routine determines if n possess a prime factor that's larger than
//  the specified primeIndex'th prime (with 2 corresponding to primeIndex 0).
//
//  For the FFTW-2, use primeIndex = 4 to get values
//  of n for which the routine is consistently efficient.
//
//  For the FFTW-3, use primeIndex = 6 to get values
//  of n for which the routine is consistently efficient.

	bool  checkLargePrimeFactors(long n)
	{
//
//  A list of the first 30 primes
//
	static long primeList[] =
	{2, 3, 5, 7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113};

	// Exit if too large of a input argument

	if(primeIndex > 30)
	{
	printf("XXXXXXXXXXX CheckLargePrimeFactor Error XXXXXXXXXXX\n");
	printf("Requsting factorization containing primes greater than 127 \n");
	printf("Use primeIndex <= 30 \n");
	printf("XXXXXXXXXXX\n");
	exit(0);
	}

    // Factor out all primes less than primeList[primeIndex]

    long primeCount;
    for(primeCount = 0; primeCount <= primeIndex; primeCount++)
    {
    while(n%primeList[primeCount] == 0){n = n/primeList[primeCount];};
    }

    if(n == 1) return false;
    return true;
}

    void setMaximalPrimeFactor(long factor)
    {
    static long primeList[] =
	{2, 3, 5, 7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113};

	primeIndex = 0;
	for(int i = 0; i < 30; i++)
	{
	if(factor  >= primeList[i]) primeIndex = i;
	}

	if(primeIndex >= 30)
	{
    printf("XXXXXXXXXXX setMaximalPrimeFactor Error XXXXXXXXXXX\n");
	printf("Requsting maximal prime factor greater than 127 \n");
	printf("or specifying a non-prime value for maximal prime factor. \n");
	printf("XXXXXXXXXXX\n");
	exit(0);
	}

    }
	long primeIndex; // Index of maximal allowed prime in factorization
};
}

#undef DEFAULT_MAXIMAL_PRIME_INDEX_
#endif

