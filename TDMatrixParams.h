#ifndef _TDMATRIX_PARAMS_H_
#define _TDMATRIX_PARAMS_H_
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
/*
 * TDMatrixParams.h - Parameters for trapdoor sampling
 * Access points:
 *
 *   TDMatrixParams prms(n, qBits, e, m, r);
 *      long n -     small dimension of A
 *      long qBits - bit-length of modulus q
 *      long e =     # of times each factor is used
 *      long m -     large dimension of A, defaults to O(n log q)
 *      long r -     samll constant, representing omega(sqrt(log m)),
 *                   default value is 4
 *
 *   The factors in this case are chosen from TDMATRIX_SMALL_FACTORS,
 *   as many as can be represneted modulus q
 *
 * An alternative form, kept only for debugging purposes is
 *
 *   TDMatrixParams(factors, m, r, qBits, q, n, e);
 */
#include <cstdio>
#include <iostream>
#include <cassert>
#include <NTL/ZZ.h>
#include <NTL/mat_lzz_p.h>

#include "mat_l.h"
#include "stash.h"


class TDMatrixParams
{
public:
    static long smallFactors[]; // A table of small co-prime factors (6-7.5 bits)
    // C++ doesn't support static initialization here, see TDMatrix.cpp
#define TDMATRIX_NUM_SMALL_FACTORS 28
#define TDMATRIX_SMALL_FACTORS                                         \
  181, 179, 173, 167, 163, 157, 151, 149, 139, 137, 131, 127,          \
  113, 109, 107, 103, 101,  97,  89,  83,  79,  73,  71, 143/*=11*13*/,\
  161/*=7*23*/, 155/*=5*31*/, 177/*=3*59*/, 134/*=2*67*/
#define TDMATRIX_QBITS_HARDWIRED 184 // product of factors is ~2^184

    long n;        // number of rows in tparamshe matrix A
    long mBar;     // =2n, # of cols in Abar, chosen at random when generating A
    long kFactors; // number of different small factors in q
    long e;        // each factor is used e times, factors^e used for CRT
    long wLen;     // = n*kFactors*e;
    long m;        // = mBar + wLen
    vec_l factors; // the small factors of q, q = \prod_i factors[i]^e
    long maxFactor;// the largest element in the factors vector
    long r;        // asymptotically \omega(srqt(log m)), practically =4

    long sigmaX;   // Gaussian parameter that we can sample with a trapdoor
    // Set to r^3 * maxFactor * sqrt(m)

    NTL::Vec<NTL::zz_pContext> zzp_context; // zz_p context for the CRT moduli
    NTL::mat_zz_p fInv;                     // f_i^{-1} (mod f_j^e) for all j>i

    // keep unused points that were sampled before, we have one stash per factor
    NTL::Vec<SampleStash> stash;


    TDMatrixParams() {} // empty constructor

    void init(long nn, long kk, long ee=3, long mm=0, long rr=4);
    TDMatrixParams(long nn, long kk, long ee=3, long mm=0, long rr=4)
    { init(nn,kk,ee,mm,rr); }

    NTL::ZZ getQ() const;

    // binary I/O, returns # of bytes read
    long writeToFile(FILE* handle) const;
    long readFromFile(FILE* handle);

    friend std::istream& operator>>(std::istream &s, TDMatrixParams& p);
};
// helper functions for setting the parameters
long getsigmaXVal(long ww, long mm, long mx, long rr=4);
long mBound(long qBits, long errBits, long sec);
long mbarBound(long qBits, long n, long sec);
long get_qBits(long kFactors, long e=3);
	   
bool operator==(const TDMatrixParams& p, const TDMatrixParams& q);
inline bool operator!=(const TDMatrixParams& p, const TDMatrixParams& q)
{
    return !(p==q);
}

// text-based I/O
std::ostream& operator<<(std::ostream &s, const TDMatrixParams& prms);
std::istream& operator>>(std::istream &s, TDMatrixParams& prms);

#endif // _TDMATRIX_PARAMS_H_
