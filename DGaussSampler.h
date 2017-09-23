#ifndef _DGAUSSSAMPLER_H_
#define _DGAUSSSAMPLER_H_
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
/***********************************************************************
 * DGaussSampler.h - Sampling from non-spherical discrete Gaussian.
 *   A vector x is sampled one entry at a time, each entry sampled from
 *   the marginal distribution, conditioned on all previous entries.
 *
 * Usage pattern:
 *
 *   Mat<long> covMat;    // initialized somehow
 *   DiscreteGaussianSampler DS(CovMat); // constructed gets covariance matrix
 *
 *   Vec<double> meanVec; // initialized somehow
 *   Vec<long> xOut;
 *   DS.SampleDiscreteGaussian(xOut,meanVec); // Sample xOut
 *
 * or:
 *
 *   DiscreteGaussianSampler DS; // empty constructor
 *   DS.InitSampler(CovMat);     // late initialization
 *   DS.SampleDiscreteGaussian(xOut,MeanVec); // sample xOut
 *************************************************************************/
#include <vector>
#include <cassert>
#include <stdexcept>
#include <limits> // defines numeric_limits<T>
#include <atomic>
#include <NTL/matrix.h>
#include <NTL/pair.h>

#include "mat_l.h"
#include "utils/timing.h"
#include "Gaussian1Dsampler.h"

#define CONTINUOUS_GAUSSIAN_SAMPLING 1

class Gaussian1Dsampler;

// A 1D Gaussian sampler, this is either a C++ wrapper around
// Martin Albrecht's C library or an alternative implementation



// Choose a small dimension-n integer matrix
void setSmall(vec_l& u, long n, double sigma);

// Choose a small n-by-m integer matrix
void setSmall(mat_l& E, long n, long m, double sigma);

/********************************************************************/
typedef NTL::Vec< NTL::Vec<double> > condSigmaVec;

// A class for sampling discrete Gaussians over Z with a convariance matrix
class DiscreteGaussianSampler
{
 private:
  int n;  //dimension of matrix
  condSigmaVec condSigmaV; // The conditional covariances
  // the first row of each conditional matrix after the next variable is chosen

  NTL::Vec< Gaussian1Dsampler* > oneDsamplers;
  // 1D samplers, one for each entry in the Gaussian vector

 public:
  // empty constructor
  DiscreteGaussianSampler() { oneDsamplers.SetLength(0); }

  // Initialize with a covariance matrix
  bool InitSampler(const NTL::Mat<long> &mCovMat); // return false on failure
  DiscreteGaussianSampler (const NTL::Mat<long> &mCovMat) {
    if (!InitSampler(mCovMat)) throw std::logic_error("Covariance not positive");
  }

  ~DiscreteGaussianSampler(); // free memory

  // Sampling routine
  void SampleDiscreteGaussian(vec_l &xOut, const NTL::Vec<double> &meanVec) const;

  long writeToFile(FILE* handle);
  long readFromFile(FILE* handle);
};

#endif // _DGAUSSSAMPLER_H_
