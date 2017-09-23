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
/*****************************************************************************
 * DGaussSampler.cpp - Sampling from non-spherical discrete Gaussian.
 *   A vector x is sampled one entry at a time, each entry sampled from
 *   the marginal distribution, conditioned on all previous entries.
 ******************************************************************************/
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
#include "Gaussian1Dsampler.h"
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

// Choose a small dimension-n integer vector
void setSmall(vec_l& u, long n, double sigma)
{
    FHE_TIMER_START;
    u.SetLength(n);
    NTL::Pair<bool,long> extra(false,0);

    // Choose a small noise vector
    Gaussian1Dsampler gSampler(sigma);
    for (long i=0; i<n; i++)
      u[i] = gSampler.getSample(/*mu=*/0.0, extra);
}


// Choose a small n-by-m integer matrix
void setSmall(mat_l& E, long n, long m, double sigma)
{
    FHE_TIMER_START;
    E.SetDims(n, m); // the noise n-by-m matrix E

    // Choose a small noise matrix
    Gaussian1Dsampler gSampler(sigma);

    EXEC_RANGE(m, first, last)

    NTL::Pair<bool,long> extra(false,0);
    for (long j=first; j<last; j++)
    for (long i=0; i<n; i++){
        //for (long j=0; j<m; j++){
	E[i][j] = gSampler.getSample(/*mu=*/0.0, extra);
      }
      EXEC_RANGE_END
}


// Initialize to a convariance matrix, returns false on failure. The i'th
// vector in condSigmaV is the top row of the conditional covariance matrix
// of entries i,...,n ocnditioned on 1,...,i-1. Also initializes a 1D sample
// with variance condSigmaV[i][i] for each entry i.


#ifdef NTL_HAVE_AVX
// ******* AVX code
// Most of this was taken verbatim from NTL's mat_lzz_p.c

#define MAT_BLK_SZ (32)


#define PAR_THRESH_SQ (200)
#define PAR_THRESH (40000)


// MUL_ADD(a, b, c): a += b*c
#ifdef NTL_HAVE_FMA
#define MUL_ADD(a, b, c) a = _mm256_fmadd_pd(b, c, a)
#else
#define MUL_ADD(a, b, c) a = _mm256_add_pd(a, _mm256_mul_pd(b, c))
#endif



static
void muladd1_by_32(double *x, const double *a, const double *b, long n)
{
   __m256d acc0=_mm256_load_pd(x + 0*4);
   __m256d acc1=_mm256_load_pd(x + 1*4);
   __m256d acc2=_mm256_load_pd(x + 2*4);
   __m256d acc3=_mm256_load_pd(x + 3*4);
   __m256d acc4=_mm256_load_pd(x + 4*4);
   __m256d acc5=_mm256_load_pd(x + 5*4);
   __m256d acc6=_mm256_load_pd(x + 6*4);
   __m256d acc7=_mm256_load_pd(x + 7*4);

   long i = 0;
   for (; i <= n-4; i +=4) {

      // the following code sequences are a bit faster than
      // just doing 4 _mm256_broadcast_sd's
      // it requires a to point to aligned storage, however

     // this one seems slightly faster
      __m256d a0101 = _mm256_broadcast_pd((const __m128d*)(a+0));
      __m256d a2323 = _mm256_broadcast_pd((const __m128d*)(a+2));


      __m256d avec0 = _mm256_permute_pd(a0101, 0);
      __m256d avec1 = _mm256_permute_pd(a0101, 0xf);
      __m256d avec2 = _mm256_permute_pd(a2323, 0);
      __m256d avec3 = _mm256_permute_pd(a2323, 0xf);

      a += 4;

      __m256d bvec;

      bvec = _mm256_load_pd(b); b += 4; MUL_ADD(acc0, avec0, bvec);
      bvec = _mm256_load_pd(b); b += 4; MUL_ADD(acc1, avec0, bvec);
      bvec = _mm256_load_pd(b); b += 4; MUL_ADD(acc2, avec0, bvec);
      bvec = _mm256_load_pd(b); b += 4; MUL_ADD(acc3, avec0, bvec);
      bvec = _mm256_load_pd(b); b += 4; MUL_ADD(acc4, avec0, bvec);
      bvec = _mm256_load_pd(b); b += 4; MUL_ADD(acc5, avec0, bvec);
      bvec = _mm256_load_pd(b); b += 4; MUL_ADD(acc6, avec0, bvec);
      bvec = _mm256_load_pd(b); b += 4; MUL_ADD(acc7, avec0, bvec);

      bvec = _mm256_load_pd(b); b += 4; MUL_ADD(acc0, avec1, bvec);
      bvec = _mm256_load_pd(b); b += 4; MUL_ADD(acc1, avec1, bvec);
      bvec = _mm256_load_pd(b); b += 4; MUL_ADD(acc2, avec1, bvec);
      bvec = _mm256_load_pd(b); b += 4; MUL_ADD(acc3, avec1, bvec);
      bvec = _mm256_load_pd(b); b += 4; MUL_ADD(acc4, avec1, bvec);
      bvec = _mm256_load_pd(b); b += 4; MUL_ADD(acc5, avec1, bvec);
      bvec = _mm256_load_pd(b); b += 4; MUL_ADD(acc6, avec1, bvec);
      bvec = _mm256_load_pd(b); b += 4; MUL_ADD(acc7, avec1, bvec);

      bvec = _mm256_load_pd(b); b += 4; MUL_ADD(acc0, avec2, bvec);
      bvec = _mm256_load_pd(b); b += 4; MUL_ADD(acc1, avec2, bvec);
      bvec = _mm256_load_pd(b); b += 4; MUL_ADD(acc2, avec2, bvec);
      bvec = _mm256_load_pd(b); b += 4; MUL_ADD(acc3, avec2, bvec);
      bvec = _mm256_load_pd(b); b += 4; MUL_ADD(acc4, avec2, bvec);
      bvec = _mm256_load_pd(b); b += 4; MUL_ADD(acc5, avec2, bvec);
      bvec = _mm256_load_pd(b); b += 4; MUL_ADD(acc6, avec2, bvec);
      bvec = _mm256_load_pd(b); b += 4; MUL_ADD(acc7, avec2, bvec);

      bvec = _mm256_load_pd(b); b += 4; MUL_ADD(acc0, avec3, bvec);
      bvec = _mm256_load_pd(b); b += 4; MUL_ADD(acc1, avec3, bvec);
      bvec = _mm256_load_pd(b); b += 4; MUL_ADD(acc2, avec3, bvec);
      bvec = _mm256_load_pd(b); b += 4; MUL_ADD(acc3, avec3, bvec);
      bvec = _mm256_load_pd(b); b += 4; MUL_ADD(acc4, avec3, bvec);
      bvec = _mm256_load_pd(b); b += 4; MUL_ADD(acc5, avec3, bvec);
      bvec = _mm256_load_pd(b); b += 4; MUL_ADD(acc6, avec3, bvec);
      bvec = _mm256_load_pd(b); b += 4; MUL_ADD(acc7, avec3, bvec);
   }

   for (; i < n; i++) {
      __m256d avec = _mm256_broadcast_sd(a); a++;
      __m256d bvec;

      bvec = _mm256_load_pd(b); b += 4; MUL_ADD(acc0, avec, bvec);
      bvec = _mm256_load_pd(b); b += 4; MUL_ADD(acc1, avec, bvec);
      bvec = _mm256_load_pd(b); b += 4; MUL_ADD(acc2, avec, bvec);
      bvec = _mm256_load_pd(b); b += 4; MUL_ADD(acc3, avec, bvec);
      bvec = _mm256_load_pd(b); b += 4; MUL_ADD(acc4, avec, bvec);
      bvec = _mm256_load_pd(b); b += 4; MUL_ADD(acc5, avec, bvec);
      bvec = _mm256_load_pd(b); b += 4; MUL_ADD(acc6, avec, bvec);
      bvec = _mm256_load_pd(b); b += 4; MUL_ADD(acc7, avec, bvec);
   }


   _mm256_store_pd(x + 0*4, acc0);
   _mm256_store_pd(x + 1*4, acc1);
   _mm256_store_pd(x + 2*4, acc2);
   _mm256_store_pd(x + 3*4, acc3);
   _mm256_store_pd(x + 4*4, acc4);
   _mm256_store_pd(x + 5*4, acc5);
   _mm256_store_pd(x + 6*4, acc6);
   _mm256_store_pd(x + 7*4, acc7);
}

// experiment: process two rows at a time
#ifndef NTL_HAVE_FMA
static
void muladd2_by_32(double *x, const double *a, const double *b, long n)
{
   __m256d avec0, avec1, bvec;
   __m256d acc00, acc01, acc02, acc03;
   __m256d acc10, acc11, acc12, acc13;

   // round 0

   acc00=_mm256_load_pd(x + 0*4 + 0*MAT_BLK_SZ);
   acc01=_mm256_load_pd(x + 1*4 + 0*MAT_BLK_SZ);
   acc02=_mm256_load_pd(x + 2*4 + 0*MAT_BLK_SZ);
   acc03=_mm256_load_pd(x + 3*4 + 0*MAT_BLK_SZ);

   acc10=_mm256_load_pd(x + 0*4 + 1*MAT_BLK_SZ);
   acc11=_mm256_load_pd(x + 1*4 + 1*MAT_BLK_SZ);
   acc12=_mm256_load_pd(x + 2*4 + 1*MAT_BLK_SZ);
   acc13=_mm256_load_pd(x + 3*4 + 1*MAT_BLK_SZ);

   for (long i = 0; i < n; i++) {
      avec0 = _mm256_broadcast_sd(&a[i]);
      avec1 = _mm256_broadcast_sd(&a[i+MAT_BLK_SZ]);

      bvec = _mm256_load_pd(&b[i*MAT_BLK_SZ+0*4]); MUL_ADD(acc00, avec0, bvec); MUL_ADD(acc10, avec1, bvec);
      bvec = _mm256_load_pd(&b[i*MAT_BLK_SZ+1*4]); MUL_ADD(acc01, avec0, bvec); MUL_ADD(acc11, avec1, bvec);
      bvec = _mm256_load_pd(&b[i*MAT_BLK_SZ+2*4]); MUL_ADD(acc02, avec0, bvec); MUL_ADD(acc12, avec1, bvec);
      bvec = _mm256_load_pd(&b[i*MAT_BLK_SZ+3*4]); MUL_ADD(acc03, avec0, bvec); MUL_ADD(acc13, avec1, bvec);
   }


   _mm256_store_pd(x + 0*4 + 0*MAT_BLK_SZ, acc00);
   _mm256_store_pd(x + 1*4 + 0*MAT_BLK_SZ, acc01);
   _mm256_store_pd(x + 2*4 + 0*MAT_BLK_SZ, acc02);
   _mm256_store_pd(x + 3*4 + 0*MAT_BLK_SZ, acc03);

   _mm256_store_pd(x + 0*4 + 1*MAT_BLK_SZ, acc10);
   _mm256_store_pd(x + 1*4 + 1*MAT_BLK_SZ, acc11);
   _mm256_store_pd(x + 2*4 + 1*MAT_BLK_SZ, acc12);
   _mm256_store_pd(x + 3*4 + 1*MAT_BLK_SZ, acc13);

   // round 1

   acc00=_mm256_load_pd(x + 4*4 + 0*MAT_BLK_SZ);
   acc01=_mm256_load_pd(x + 5*4 + 0*MAT_BLK_SZ);
   acc02=_mm256_load_pd(x + 6*4 + 0*MAT_BLK_SZ);
   acc03=_mm256_load_pd(x + 7*4 + 0*MAT_BLK_SZ);

   acc10=_mm256_load_pd(x + 4*4 + 1*MAT_BLK_SZ);
   acc11=_mm256_load_pd(x + 5*4 + 1*MAT_BLK_SZ);
   acc12=_mm256_load_pd(x + 6*4 + 1*MAT_BLK_SZ);
   acc13=_mm256_load_pd(x + 7*4 + 1*MAT_BLK_SZ);

   for (long i = 0; i < n; i++) {
      avec0 = _mm256_broadcast_sd(&a[i]);
      avec1 = _mm256_broadcast_sd(&a[i+MAT_BLK_SZ]);

      bvec = _mm256_load_pd(&b[i*MAT_BLK_SZ+0*4+MAT_BLK_SZ/2]); MUL_ADD(acc00, avec0, bvec); MUL_ADD(acc10, avec1, bvec);
      bvec = _mm256_load_pd(&b[i*MAT_BLK_SZ+1*4+MAT_BLK_SZ/2]); MUL_ADD(acc01, avec0, bvec); MUL_ADD(acc11, avec1, bvec);
      bvec = _mm256_load_pd(&b[i*MAT_BLK_SZ+2*4+MAT_BLK_SZ/2]); MUL_ADD(acc02, avec0, bvec); MUL_ADD(acc12, avec1, bvec);
      bvec = _mm256_load_pd(&b[i*MAT_BLK_SZ+3*4+MAT_BLK_SZ/2]); MUL_ADD(acc03, avec0, bvec); MUL_ADD(acc13, avec1, bvec);
   }


   _mm256_store_pd(x + 4*4 + 0*MAT_BLK_SZ, acc00);
   _mm256_store_pd(x + 5*4 + 0*MAT_BLK_SZ, acc01);
   _mm256_store_pd(x + 6*4 + 0*MAT_BLK_SZ, acc02);
   _mm256_store_pd(x + 7*4 + 0*MAT_BLK_SZ, acc03);

   _mm256_store_pd(x + 4*4 + 1*MAT_BLK_SZ, acc10);
   _mm256_store_pd(x + 5*4 + 1*MAT_BLK_SZ, acc11);
   _mm256_store_pd(x + 6*4 + 1*MAT_BLK_SZ, acc12);
   _mm256_store_pd(x + 7*4 + 1*MAT_BLK_SZ, acc13);

}
#endif

// experiment: process three rows at a time
// NOTE: this makes things slower on an AVX1 platform --- not enough registers
// it could be faster on AVX2/FMA, where there should be enough registers

static
void muladd3_by_32(double *x, const double *a, const double *b, long n)
{
   __m256d avec0, avec1, avec2, bvec;
   __m256d acc00, acc01, acc02, acc03;
   __m256d acc10, acc11, acc12, acc13;
   __m256d acc20, acc21, acc22, acc23;


   // round 0

   acc00=_mm256_load_pd(x + 0*4 + 0*MAT_BLK_SZ);
   acc01=_mm256_load_pd(x + 1*4 + 0*MAT_BLK_SZ);
   acc02=_mm256_load_pd(x + 2*4 + 0*MAT_BLK_SZ);
   acc03=_mm256_load_pd(x + 3*4 + 0*MAT_BLK_SZ);

   acc10=_mm256_load_pd(x + 0*4 + 1*MAT_BLK_SZ);
   acc11=_mm256_load_pd(x + 1*4 + 1*MAT_BLK_SZ);
   acc12=_mm256_load_pd(x + 2*4 + 1*MAT_BLK_SZ);
   acc13=_mm256_load_pd(x + 3*4 + 1*MAT_BLK_SZ);

   acc20=_mm256_load_pd(x + 0*4 + 2*MAT_BLK_SZ);
   acc21=_mm256_load_pd(x + 1*4 + 2*MAT_BLK_SZ);
   acc22=_mm256_load_pd(x + 2*4 + 2*MAT_BLK_SZ);
   acc23=_mm256_load_pd(x + 3*4 + 2*MAT_BLK_SZ);

   for (long i = 0; i < n; i++) {
      avec0 = _mm256_broadcast_sd(&a[i]);
      avec1 = _mm256_broadcast_sd(&a[i+MAT_BLK_SZ]);
      avec2 = _mm256_broadcast_sd(&a[i+2*MAT_BLK_SZ]);

      bvec = _mm256_load_pd(&b[i*MAT_BLK_SZ+0*4]); MUL_ADD(acc00, avec0, bvec); MUL_ADD(acc10, avec1, bvec); MUL_ADD(acc20, avec2, bvec);
      bvec = _mm256_load_pd(&b[i*MAT_BLK_SZ+1*4]); MUL_ADD(acc01, avec0, bvec); MUL_ADD(acc11, avec1, bvec); MUL_ADD(acc21, avec2, bvec);
      bvec = _mm256_load_pd(&b[i*MAT_BLK_SZ+2*4]); MUL_ADD(acc02, avec0, bvec); MUL_ADD(acc12, avec1, bvec); MUL_ADD(acc22, avec2, bvec);
      bvec = _mm256_load_pd(&b[i*MAT_BLK_SZ+3*4]); MUL_ADD(acc03, avec0, bvec); MUL_ADD(acc13, avec1, bvec); MUL_ADD(acc23, avec2, bvec);
   }


   _mm256_store_pd(x + 0*4 + 0*MAT_BLK_SZ, acc00);
   _mm256_store_pd(x + 1*4 + 0*MAT_BLK_SZ, acc01);
   _mm256_store_pd(x + 2*4 + 0*MAT_BLK_SZ, acc02);
   _mm256_store_pd(x + 3*4 + 0*MAT_BLK_SZ, acc03);

   _mm256_store_pd(x + 0*4 + 1*MAT_BLK_SZ, acc10);
   _mm256_store_pd(x + 1*4 + 1*MAT_BLK_SZ, acc11);
   _mm256_store_pd(x + 2*4 + 1*MAT_BLK_SZ, acc12);
   _mm256_store_pd(x + 3*4 + 1*MAT_BLK_SZ, acc13);

   _mm256_store_pd(x + 0*4 + 2*MAT_BLK_SZ, acc20);
   _mm256_store_pd(x + 1*4 + 2*MAT_BLK_SZ, acc21);
   _mm256_store_pd(x + 2*4 + 2*MAT_BLK_SZ, acc22);
   _mm256_store_pd(x + 3*4 + 2*MAT_BLK_SZ, acc23);

   // round 1

   acc00=_mm256_load_pd(x + 4*4 + 0*MAT_BLK_SZ);
   acc01=_mm256_load_pd(x + 5*4 + 0*MAT_BLK_SZ);
   acc02=_mm256_load_pd(x + 6*4 + 0*MAT_BLK_SZ);
   acc03=_mm256_load_pd(x + 7*4 + 0*MAT_BLK_SZ);

   acc10=_mm256_load_pd(x + 4*4 + 1*MAT_BLK_SZ);
   acc11=_mm256_load_pd(x + 5*4 + 1*MAT_BLK_SZ);
   acc12=_mm256_load_pd(x + 6*4 + 1*MAT_BLK_SZ);
   acc13=_mm256_load_pd(x + 7*4 + 1*MAT_BLK_SZ);

   acc20=_mm256_load_pd(x + 4*4 + 2*MAT_BLK_SZ);
   acc21=_mm256_load_pd(x + 5*4 + 2*MAT_BLK_SZ);
   acc22=_mm256_load_pd(x + 6*4 + 2*MAT_BLK_SZ);
   acc23=_mm256_load_pd(x + 7*4 + 2*MAT_BLK_SZ);

   for (long i = 0; i < n; i++) {
      avec0 = _mm256_broadcast_sd(&a[i]);
      avec1 = _mm256_broadcast_sd(&a[i+MAT_BLK_SZ]);
      avec2 = _mm256_broadcast_sd(&a[i+2*MAT_BLK_SZ]);

      bvec = _mm256_load_pd(&b[i*MAT_BLK_SZ+0*4+MAT_BLK_SZ/2]); MUL_ADD(acc00, avec0, bvec); MUL_ADD(acc10, avec1, bvec); MUL_ADD(acc20, avec2, bvec);
      bvec = _mm256_load_pd(&b[i*MAT_BLK_SZ+1*4+MAT_BLK_SZ/2]); MUL_ADD(acc01, avec0, bvec); MUL_ADD(acc11, avec1, bvec); MUL_ADD(acc21, avec2, bvec);
      bvec = _mm256_load_pd(&b[i*MAT_BLK_SZ+2*4+MAT_BLK_SZ/2]); MUL_ADD(acc02, avec0, bvec); MUL_ADD(acc12, avec1, bvec); MUL_ADD(acc22, avec2, bvec);
      bvec = _mm256_load_pd(&b[i*MAT_BLK_SZ+3*4+MAT_BLK_SZ/2]); MUL_ADD(acc03, avec0, bvec); MUL_ADD(acc13, avec1, bvec); MUL_ADD(acc23, avec2, bvec);
   }


   _mm256_store_pd(x + 4*4 + 0*MAT_BLK_SZ, acc00);
   _mm256_store_pd(x + 5*4 + 0*MAT_BLK_SZ, acc01);
   _mm256_store_pd(x + 6*4 + 0*MAT_BLK_SZ, acc02);
   _mm256_store_pd(x + 7*4 + 0*MAT_BLK_SZ, acc03);

   _mm256_store_pd(x + 4*4 + 1*MAT_BLK_SZ, acc10);
   _mm256_store_pd(x + 5*4 + 1*MAT_BLK_SZ, acc11);
   _mm256_store_pd(x + 6*4 + 1*MAT_BLK_SZ, acc12);
   _mm256_store_pd(x + 7*4 + 1*MAT_BLK_SZ, acc13);

   _mm256_store_pd(x + 4*4 + 2*MAT_BLK_SZ, acc20);
   _mm256_store_pd(x + 5*4 + 2*MAT_BLK_SZ, acc21);
   _mm256_store_pd(x + 6*4 + 2*MAT_BLK_SZ, acc22);
   _mm256_store_pd(x + 7*4 + 2*MAT_BLK_SZ, acc23);

}


static inline
void muladd_all_by_32(long first, long last, double *x, const double *a, const double *b, long n)
{
   long i = first;
#ifdef NTL_HAVE_FMA
   // processing three rows at a time is faster
   for (; i <= last-3; i+=3)
      muladd3_by_32(x + i*MAT_BLK_SZ, a + i*MAT_BLK_SZ, b, n);
   for (; i < last; i++)
      muladd1_by_32(x + i*MAT_BLK_SZ, a + i*MAT_BLK_SZ, b, n);
#else
   // process only two rows at a time: not enough registers :-(
   for (; i <= last-2; i+=2)
      muladd2_by_32(x + i*MAT_BLK_SZ, a + i*MAT_BLK_SZ, b, n);
   for (; i < last; i++)
      muladd1_by_32(x + i*MAT_BLK_SZ, a + i*MAT_BLK_SZ, b, n);
#endif
}

// this assumes n is a multiple of 16
static inline
void muladd_interval(double * NTL_RESTRICT x, double * NTL_RESTRICT y, double c, long n)
{
   __m256d xvec0, xvec1, xvec2, xvec3;
   __m256d yvec0, yvec1, yvec2, yvec3;

   __m256d cvec = _mm256_broadcast_sd(&c);

   for (long i = 0; i < n; i += 16, x += 16, y += 16) {
      xvec0 = _mm256_load_pd(x+0*4);
      xvec1 = _mm256_load_pd(x+1*4);
      xvec2 = _mm256_load_pd(x+2*4);
      xvec3 = _mm256_load_pd(x+3*4);

      yvec0 = _mm256_load_pd(y+0*4);
      yvec1 = _mm256_load_pd(y+1*4);
      yvec2 = _mm256_load_pd(y+2*4);
      yvec3 = _mm256_load_pd(y+3*4);

      MUL_ADD(xvec0, yvec0, cvec);
      MUL_ADD(xvec1, yvec1, cvec);
      MUL_ADD(xvec2, yvec2, cvec);
      MUL_ADD(xvec3, yvec3, cvec);

      _mm256_store_pd(x + 0*4, xvec0);
      _mm256_store_pd(x + 1*4, xvec1);
      _mm256_store_pd(x + 2*4, xvec2);
      _mm256_store_pd(x + 3*4, xvec3);
   }
}


static
bool Alt_sampler(Vec< Vec<double> >& X, const Mat<long>& A)
{
   FHE_TIMER_START;

   long n = A.NumRows();

   if (A.NumCols() != n)
      LogicError("sampler: nonsquare matrix");

   if (NTL_OVERFLOW(n, MAT_BLK_SZ, 0)) ResourceError("dimension too large");

   long npanels = (n+MAT_BLK_SZ-1)/MAT_BLK_SZ;

   Vec< AlignedArray<double> > M;
   M.SetLength(npanels);
   for (long panel = 0; panel < npanels; panel++) {
      M[panel].SetLength(n*MAT_BLK_SZ);
      double *panelp = &M[panel][0];

      for (long r = 0; r < n*MAT_BLK_SZ; r++) panelp[r] = 0;
   }

   // copy A into panels
   for (long jj = 0, panel = 0; jj < n; jj += MAT_BLK_SZ, panel++) {
      long j_max = min(jj+MAT_BLK_SZ, n);
      double *panelp = &M[panel][0];

      for (long i = 0; i < n; i++, panelp += MAT_BLK_SZ) {
         const long *ap = A[i].elts() + jj;

         for (long j = jj; j < j_max; j++)
            panelp[j-jj] = ap[j-jj];
      }
   }

   X.SetLength(n);

   for (long kk = 0, kpanel = 0; kk < n; kk += MAT_BLK_SZ, kpanel++) {
      long k_max = min(kk+MAT_BLK_SZ, n);

      double * NTL_RESTRICT kpanelp = &M[kpanel][0];

      for (long k = kk; k < k_max; k++) {

         X[k].SetLength(n-k);
         double * NTL_RESTRICT Xk = X[k].elts();
         for (long i = k; i < n; i++)
            Xk[i-k] = kpanelp[i*MAT_BLK_SZ + (k-kk)];

         double * NTL_RESTRICT y = &kpanelp[k*MAT_BLK_SZ];
         double pivot = y[k-kk];
         y[k-kk] = 1;

         if (pivot <= 0) {
            return false;
         }

         double pivot_inv = 1/pivot;

         for (long i = k+1; i < n; i++) {
            double * NTL_RESTRICT x = &kpanelp[i*MAT_BLK_SZ];
            double t1 = -x[k-kk]*pivot_inv;
            x[k-kk] = 0;
            muladd_interval(x, y, t1, MAT_BLK_SZ);
         }

         for (long j = k; j < k_max; j++) y[j-kk] = 0;
      }


      // finished processing current kpanel
      // next, reduce and apply to all other kpanels

      bool seq = double(npanels-(kpanel+1))*double(n)*double(MAT_BLK_SZ)*double(MAT_BLK_SZ) < PAR_THRESH;

      NTL_GEXEC_RANGE(seq, npanels-(kpanel+1), first, last)
      NTL_IMPORT(n)
      NTL_IMPORT(kpanel)
      NTL_IMPORT(kpanelp)
      NTL_IMPORT(kk)
      NTL_IMPORT(k_max)

      AlignedArray<double> buf_store;
      buf_store.SetLength(MAT_BLK_SZ*MAT_BLK_SZ);
      double *buf = &buf_store[0];

      for (long index = first; index < last; index++) {
         long jpanel = index + kpanel+1;

         double * NTL_RESTRICT jpanelp = &M[jpanel][0];

         // copy block number kpanel (the one on the diagonal)  into buf

         for (long i = 0; i < (k_max-kk)*MAT_BLK_SZ; i++)
            buf[i] = jpanelp[kk*MAT_BLK_SZ+i];

         // jpanel += kpanel*buf

         muladd_all_by_32(kk, n, jpanelp, kpanelp, buf, k_max-kk);
      }

      NTL_GEXEC_RANGE_END

   }

   return true;

}
#endif

//create the nested covariance matrices
static
bool Plain_sampler(Vec< Vec<double> >& X, const Mat<long>& A)
{
    FHE_TIMER_START;

    if (A[0][0] < 0)   //sigma too small
    {
        return false; // failed, variance cannot be negative
    }

    long n = A.NumRows();

    if (A.NumCols() != n)
       LogicError("sampler: nonsquare matrix");

    X.SetLength(n); // allocate space for n vectors

    Mat<double> M;
    conv(M, A); // convert to floating point

    // 1st vector of the conditional covariance = 1st column of covariance
    X[0].SetLength(n);
    for (long k = 0; k < n; k++)
        X[0][k] = M[k][0];

    // Compute the rest of the conditional covariance, one vector at a time


    for (long k = 1; k<n; k++)
    {
        // create the new covariance: sig y|x = sig_yy - sig_yx *sig_xx^{-1} *sig_xy


        //get the subset covariance matrix

        NTL_GEXEC_RANGE(n-k < 100, n-k, first, last)

        NTL_IMPORT(n)
        NTL_IMPORT(k)

        double * NTL_RESTRICT M_0 = &M[k-1][0];
        double pivot = M_0[k-1];
        double pivot_inv = 1.0/pivot;

        for (long row = first; row < last; row++) {
            double * NTL_RESTRICT M_r = &M[k+row][0];
            double fac = -M_r[k-1]*pivot_inv;

            for (long col = k; col < n; col++) {
                M_r[col] += fac * (M_0[col]);
            }
        }

        NTL_GEXEC_RANGE_END

        if (M[k][k] < 0)
        {
            return false; //failed, sigma too small
        }

        // copy the result to the next conditional covariance vector
        X[k].SetLength(n-k);
        for (long i = 0; i < n-k; i++)
            X[k][i] = M[k+i][k];
    }

    return true;
}



bool DiscreteGaussianSampler::InitSampler(const Mat<long> &covMat)
{
    FHE_TIMER_START;

    n = covMat.NumRows();

    bool res = false;

#ifdef NTL_HAVE_AVX
#warning "using AVX code"
   if (n >= 512) {
      res = Alt_sampler(condSigmaV, covMat);
   }
   else
#endif
   {
      res = Plain_sampler(condSigmaV, covMat);
   }

   if (!res) return false;

    //clear if there are previous 1D-samplers
    for (long i = 0; i < oneDsamplers.length(); i++)
        if (oneDsamplers[i] != NULL) delete oneDsamplers[i];

    // create table of 1D-samplers
    oneDsamplers.SetLength(n);

    for (long iVar=0; iVar<n; iVar++)
    {
        if (condSigmaV[iVar][0] < 0)
        {
            cout << "conditional sigma too small" << endl;
            return false; //conditional sigma too small
            }

        oneDsamplers[iVar] =            //create a 1D sampler for this sigma
        (Gaussian1Dsampler*) new Gaussian1Dsampler(sqrt(condSigmaV[iVar][0]));

    }

    return true;
}



// FIXME: better to use UniquePtr or something to manage the entries
// in oneDsamplers, rather than new/delete


//delete the array of sampler instances
DiscreteGaussianSampler::~DiscreteGaussianSampler()
{
    FHE_TIMER_START;
    for (long i = 0; i < oneDsamplers.length(); i++)
        if (oneDsamplers[i]) delete oneDsamplers[i];
}

// Discrete sampling according to the conditional variance and mean. Choose
//each entry in x depending on the conditioned distribution of previous ones.

void DiscreteGaussianSampler::SampleDiscreteGaussian(vec_l& xOut, const Vec<double>& meanVec) const
{
    FHE_TIMER_START;

    //for this val, find the rest of the values

    //x,y are matrices
    //m y|x = m_y + sig_yx * sig ^-1_xx (x - m_x)
    //multiple matrices

    // scratch variables for internal computation
    Vec<double> marginalMeanVec = meanVec;
    Vec<double> nextMarginalMeanVec;

    // work with pointers for easier swapping
    Vec<double>* curMeanVec = &marginalMeanVec;
    Vec<double>* nextMeanVec = &nextMarginalMeanVec;

    //find marginal distribution for each variable

    xOut.SetLength(n);
    for (long iVar = 0; iVar<n; iVar++)
    {
        //get the next sample
        double mean = (*curMeanVec)[0];
        xOut[iVar] = oneDsamplers[iVar]->getSample(mean);

        // update the mean vector:
        //    m(i|j) = m(i) + sigma(i,j)*sigma(j,j)^-1*(x(j) - m(j))
        // only need the first element for the next item

        // Next mean vector has one dimension less
        nextMeanVec->SetLength(curMeanVec->length() -1);

        double diff = xOut[iVar] - mean;  // deviation from mean

        // Look at the next vector from the conditional covariance
        const Vec<double>& covVec = condSigmaV[iVar];

        // newMean = curMean[1..n-1] + covVec * diff / sigma_XX

        double invSig = 1.0 / covVec[0]; // 1/sigma_XX
        for (long row = 0; row<nextMeanVec->length(); row++)
            (*nextMeanVec)[row]
                = (*curMeanVec)[row+1] + (covVec[row+1] *invSig) *diff;

        // Set the updated mean vector as the current one
        swap(nextMeanVec,curMeanVec); // pointer swap
    }
}



long DiscreteGaussianSampler::writeToFile(FILE* handle)
{
    FHE_TIMER_START;
    long count = fwrite(&n,sizeof(long),1,handle); // write dimension

    long nConvSig = condSigmaV.length();
    count += fwrite(&nConvSig,sizeof(long),1,handle);

    for (long i = 0; i < nConvSig; i++)
    {
    long nconvSigmaVi = condSigmaV[i].length();
    count += fwrite(&nconvSigmaVi,sizeof(long),1,handle);
    count += fwrite(condSigmaV[i].elts(),condSigmaV[i].length()*sizeof(double),1,handle);
    }

    long oLength = oneDsamplers.length();
    count += fwrite(&oLength,sizeof(long),1,handle);


    for (long i = 0; i < oLength; i++)
        count+= oneDsamplers[i]->writeToFile(handle);

    return count;
}

long DiscreteGaussianSampler::readFromFile(FILE* handle)
{
    FHE_TIMER_START;
    long count = fread(&n,sizeof(long),1,handle); // read dimension

    long nConvSig;
    count += fread(&nConvSig,sizeof(long),1,handle);

    condSigmaV.SetLength(nConvSig);

    for (long i = 0; i < nConvSig; i++)
    {
    long nconvSigmaVi;
    count += fread(&nconvSigmaVi,sizeof(long),1,handle);
    condSigmaV[i].SetLength(nconvSigmaVi);
    count += fread(condSigmaV[i].elts(),nconvSigmaVi*sizeof(double),1,handle);
    }

    long oLength;
    count += fread(&oLength,sizeof(long),1,handle);
    oneDsamplers.SetLength(oLength);

    for (long iVar=0; iVar<oLength; iVar++)
    {
        oneDsamplers[iVar] =            //create a 1D sampler for this sigma
            (Gaussian1Dsampler*) new Gaussian1Dsampler(sqrt(condSigmaV[iVar][0]));
    }

    for (long i = 0; i < oLength; i++)
        count+= oneDsamplers[i]->readFromFile(handle);

    return count;
}
