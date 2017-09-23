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
#include <NTL/ZZ.h>
#include <NTL/mat_lzz_p.h>
#include <limits>
#include <NTL/BasicThreadPool.h>
#include <NTL/ZZ.h>
#include <NTL/FFT.h>
#include <NTL/SmartPtr.h>

NTL_CLIENT
#include "mat_l.h"
#include "vec_l.h"
#include "DGaussSampler.h"
#include "utils/tools.h"


#ifdef NTL_HAVE_AVX
#warning "HAVE_AVX"

#include <immintrin.h>

#ifdef NTL_HAVE_FMA
#define MUL_ADD(a, b, c) a = _mm256_fmadd_pd(b, c, a)
#else
#define MUL_ADD(a, b, c) a = _mm256_add_pd(a, _mm256_mul_pd(b, c))
#endif
#endif


// Two static variables for debugging purposes
std::atomic<double> Gaussian1Dsampler::maxSigma(0.0);// keep largest stdev we've seen
std::atomic<int> Gaussian1Dsampler::maxSample(0);   // keep largest sample ever drawn


long Gaussian1Dsampler::writeToFile(FILE* handle)
{
    FHE_TIMER_START;
    long count = fwrite(&sigma,sizeof(sigma),1,handle);

    return count;
}

long Gaussian1Dsampler::readFromFile(FILE* handle)
{
    FHE_TIMER_START;
    long count = fread(&sigma,sizeof(sigma),1,handle);

    return count;
}

#ifdef CONTINUOUS_GAUSSIAN_SAMPLING
/* generateGaussianNoiseMars uses Marsaglia polar method
 * to generate Gaussian samples
 */
NTL::Pair<double,double>
Gaussian1Dsampler::generateGaussianNoiseMars(double mean, double stdDev)
{
    FHE_TIMER_START;
    double u, v, s;
    do    // choose (u,v) in the square [-1,1] x [-1,1]
    {
        u = (NTL::RandomBnd(NTL_SP_BOUND)/((double)NTL_SP_BOUND)) *2.0 -1.0;
        v = (NTL::RandomBnd(NTL_SP_BOUND)/((double)NTL_SP_BOUND)) *2.0 -1.0;
        s = u * u + v * v;
    }
    while( (s >= 1.0) || (s == 0.0) );// until you find a pair in the unit sphere

    s = sqrt(-2.0 * log(s) / s);
    double sample1 = mean + stdDev * u * s;
    double sample2 = mean + stdDev * v * s;

    return NTL::Pair<double,double>(sample1,sample2);
}

/* generateGaussianNoiseBox uses Box-Muller to generate Gaussian Samples
 */
NTL::Pair<double,double>
Gaussian1Dsampler::generateGaussianNoiseBox(double mu, double stDev)
{
    FHE_TIMER_START;
    const double epsilon = std::numeric_limits<double>::min();
    const double two_pi = 2.0*3.14159265358979323846;

    double u1, u2;
    do
    {
        u1 = ((NTL::RandomBnd(NTL_SP_BOUND))/((double)NTL_SP_BOUND));
        u2 = ((NTL::RandomBnd(NTL_SP_BOUND))/((double)NTL_SP_BOUND));
    }
    while ( u1 <= epsilon );

    double z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
    double z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);

    double sample1 = z0 * stDev + mu;
    double sample2 = z1 * stDev + mu;

    return NTL::Pair<double,double>(sample1,sample2);
}
#endif // ifdef CONTINUOUS_GAUSSIAN_SAMPLING
