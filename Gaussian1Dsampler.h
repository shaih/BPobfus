#ifndef _GAUSSIAN1DSAMPLER_H_
#define _GAUSSIAN1DSAMPLER_H_
/* Copyright (C) 2017 IBM Corp.
 *  Licensed under the Apache License, Version 2.0 (the "License"); 
 * you may not use this file except in compliance with the License. 
 * You may obtain a copy of the License at
 *     http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, 
 * software distributed under the License is distributed on an
 * "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
 * either express or implied. See the License for the specific
 * language governing permissions and limitations under the License. 
 */
/*******************************************************************
Gaussian1Dsampler: defines a 1D pseudo-random number sampler, either
using Martin Albresht's online library or implementing two optional
methods, Box–Muller transform Marsaglia polar method
********************************************************************/
#include <vector>
#include <cassert>
#include <stdexcept>
#include <limits> // defines numeric_limits<T>
#include <atomic>
#include <NTL/matrix.h>
#include <NTL/pair.h>

#include "mat_l.h"
#include "utils/timing.h"

// Define CONTINUOUS_GAUSSIAN_SAMPLING 1 to get Marsaglia polar method.
//        CONTINUOUS_GAUSSIAN_SAMPLING 2 to get Box-Muller method,
// Keep it undefined to get Martin Albresht's library
#define CONTINUOUS_GAUSSIAN_SAMPLING 1


// A 1D Gaussian sampler, this is either a C++ wrapper around
// Martin Albrecht's C library or an alternative implementation


#ifndef CONTINUOUS_GAUSSIAN_SAMPLING // Martin's library
#include <mpfr.h>

class Gaussian1Dsampler
{
  double sigma;
  dgs_disc_gauss_dp_t* self;

 public:
  // Two static variables for debugging purposes
  static std::atomic<double> maxSigma; // keep the largest stdev we've seen
  static std::atomic<int> maxSample;  // keep the largest sample ever drawn

  Gaussian1Dsampler(double sig, long c=0, long tau=6,
	   dgs_disc_gauss_alg_t alg=DGS_DISC_GAUSS_UNIFORM_ONLINE)
    {
      assert(sig >= 0.0);
      sigma = sig;
      if (sig>maxSigma) maxSigma.store(sig); // record largest sigma for debugging
      self = dgs_disc_gauss_dp_init(sig, c, tau, alg);
    }

  ~Gaussian1Dsampler() { if (self) dgs_disc_gauss_dp_clear(self); }

  long getSample(double mu=0.0) {
    FHE_TIMER_START;
    long ret = self->call(self);
    if (ret>maxSample) maxSample.store(ret);
    return ret;
  }
};

#else  // ifdef CONTINUOUS_GAUSSIAN_SAMPLING
class Gaussian1Dsampler
{
  double sigma;

  // Two different methods for generating continuous normal R.V.
  // Both methods return two such normal R.V.'s on each call
  static NTL::Pair<double,double>
    generateGaussianNoiseMars(double mu,  double sigma);
  static NTL::Pair<double,double>
    generateGaussianNoiseBox(double mu, double sigma);

 public:
  // Two static variables for debugging purposes
  static std::atomic<double> maxSigma; // keep the largest stdev we've seen
  static std::atomic<int> maxSample;  // keep the largest sample ever drawn

  Gaussian1Dsampler(double sig) {
    assert(sig >= 0.0);
    sigma = sig;

    if (sig>maxSigma) maxSigma.store(sig); // record largest sigma for debugging
  }

  long getSample(double mu=0.0) {
    FHE_TIMER_START;
#if CONTINUOUS_GAUSSIAN_SAMPLING==1
    NTL::Pair<double,double> samples = generateGaussianNoiseMars(mu,sigma);
#else
    NTL::Pair<double,double> samples = generateGaussianNoiseBox(mu,sigma);
#endif
    return (long)round(samples.a);
  }
  // A convenience method that keeps 2nd sample in the extra argument
  long getSample(double mu, NTL::Pair<bool,long>& extra) {
    FHE_TIMER_START;
    if (extra.a) { extra.a=false; return extra.b; }
    else { 
#if CONTINUOUS_GAUSSIAN_SAMPLING==1
      NTL::Pair<double,double> samples = generateGaussianNoiseMars(mu,sigma);
#else
      NTL::Pair<double,double> samples = generateGaussianNoiseBox(mu,sigma);
#endif
      extra.a = true; extra.b = (long)round(samples.b);
      return (long)round(samples.a);
    }
  }

  long writeToFile(FILE* handle);
  long readFromFile(FILE* handle);
};
#endif // CONTINUOUS_GAUSSIAN_SAMPLING

#endif // _GAUSSIAN1DSAMPLER_H_
