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

/* TDMatrixParams.cpp - Parameters for trapdoor sampling
 */
#include <stdexcept>
#include <NTL/vector.h>
#include "utils/tools.h"
#include "utils/timing.h"
#include "TDMatrixParams.h"
#include <stdexcept>
#include <sys/stat.h>
#include <unistd.h>

//#define DEBUGPRINT
//#define DEBUG

NTL_CLIENT

// A static table of small co-prime factors (6-7.5 bits)
long TDMatrixParams::smallFactors[TDMATRIX_NUM_SMALL_FACTORS]
    = { TDMATRIX_SMALL_FACTORS }; // TDMATRIX_SMALL_FACTORS list in TDMatrix.h


// Calculate the sigmaX value. This value must be large enough so
// that sigmaX*I - r*maxFactor * (R/I)*(R^t|I) is positive definite.
// We compute an upper bound s on the singular values of R, then set
// sigmaX = r*maxFactor * s^2.
long getsigmaXVal(long wLen, long mBar, long maxFactor, long r)
{
  double s = r*(sqrt(wLen)+sqrt(mBar)+6);
  double sigmaZ = r*maxFactor;
  long val = ceil(sigmaZ*sigmaZ * s*s);

#ifdef DEBUG
  cout << "sigmaX val=" << val << endl;
#endif // DEBUG

  return val;
}

// lower bound the dimension to get sec bits of security
long mBound(long qBits, long errBits, long sec)
{
  return ceil((qBits-errBits)*(sec+110)/7.2);
}

long mbarBound(long qBits, long n, long sec)
{
  return ceil((2+sqrt(sec))*sqrt(n*qBits));
}

long get_qBits(long kFactors, long e)
{
  assert(kFactors>0 && e>0);
  if (kFactors > TDMATRIX_NUM_SMALL_FACTORS)
    kFactors = TDMATRIX_NUM_SMALL_FACTORS;

  double logProd = log(TDMatrixParams::smallFactors[0]);
  for (long i=1; i<kFactors; i++)
    logProd += log(TDMatrixParams::smallFactors[i]);

  return ceil(logProd*e/log(2.0));
}

// For each factor f_i, initialize zz_p context for F_i=f_i^e, and also
// compute f_i^{-1} mod F_j for all i<j. Returns max the largest factor
static long setFactors(mat_zz_p& fInv, Vec<zz_pContext>& zzp_context,
                       const vec_l factors, long e)
{
    FHE_TIMER_START;
    long kFactors = factors.length();

    // For each factor f_i, initialize zz_p context for F_i=f_i^e,
    // and also compute f_i^{-1} mod F_j for all i<j
    zzp_context.SetLength(kFactors);
    fInv.SetDims(kFactors,kFactors);

    NTL::zz_pPush ppush; // backup NTL's current modulus
    long maxFactor = 0;
    for (long j=kFactors-1; j>=0; j--)
    {
        if (factors[j]>maxFactor) maxFactor = factors[j];

        long f2e = NTL::power_long(factors[j],e);
        // FIXME: check that this fits in a single-precision integer
        NTL::zz_p::init(f2e);
        zzp_context[j].save(); // save the current zz_p::modulus()
        for (long i=j-1; i>=0; i--)
        {
            fInv[i][j] = NTL::inv(conv<zz_p>(factors[i])); // f_i^{-1} mod F_j
        }
    }
    return maxFactor;
}

// Initialize the parameters for a modulus q with kk factors,
// each factor is an e-th power of an even smaller number.
void TDMatrixParams::init(long nn, long kk, long ee, long mm, long rr)
{
    FHE_TIMER_START;

    assert(kk <= TDMATRIX_NUM_SMALL_FACTORS);
    factors.SetLength(kk);
    for (long i=0; i<kk; i++)
      factors[i] = smallFactors[i];
    kFactors = kk;

    e = ee;
    r = rr;
    n = nn;
    wLen = nn*kk*ee;

    if (mm > 0) { // caller supplied a value for m
      m = mm;
      mBar = m - wLen;      
      if (mBar < nn) { // check that mBar is not too tiny
        mBar = nn;
	m = mBar + wLen;
      }
    }
    else { // compute m, mBar using security formulas
      long qBits = get_qBits(kk,ee);
      m = mBound(qBits, /*errBits=*/7, /*sec=*/80);
      mBar = mbarBound(qBits, n, /*sec=*/80);
      if (m < mBar+wLen) m = mBar+wLen;
      else               mBar = m - wLen;
    }

    // For each factor f_i, initialize zz_p context for F_i=f_i^e,
    // and also compute f_i^{-1} mod F_j for all i<j
    maxFactor = setFactors(fInv, zzp_context, factors, e);

    sigmaX = getsigmaXVal(wLen, mBar, maxFactor, r);//ceil(2*r*r*r*m*maxFactor);
    // The Gaussian parameter from which we can sample with a trapdoor

    //initialize stash
    stash.SetLength(kFactors);
    for (long i = 0; i < kFactors; i++)
      stash[i].init(factors[i]);

    cout << "TDMatrixParams::init: sigmaX="<<sigmaX<<", mBar="<<mBar;
    cout << ", |q|="<<NTL::NumBits(this->getQ())<<endl;
}

//returns the produce of all the factors to the power of e, output = product(factors)^e
ZZ TDMatrixParams::getQ() const
{
    FHE_TIMER_START;
    ZZ q = to_ZZ(1L);
    for (long i=0; i<factors.length(); i++) // product of factors
        q *= factors[i];

    return NTL::power(q,e);  // return product^e
}

//outputs the different parameters in the class p into the stream s
ostream& operator<<(ostream &s, const TDMatrixParams& p)
{
    s << "[" << p.n        << " "
      << p.m        << " "
      << p.mBar     << " "
      << p.e        << " "
      << p.r        << "\n "
      << p.factors  << "]";
    return s;
}

//This function gets an input streams and sets the different parameters in p from this stream
istream& operator>>(istream &s, TDMatrixParams& p)
{
    seekPastChar(s, '[');  // this function is defined in tools.cpp
    s >> p.n;
    s >> p.m;
    s >> p.mBar;
    s >> p.e;
    s >> p.r;
    s >> p.factors;

    p.wLen = p.n * p.kFactors * p.e;

    p.maxFactor = setFactors(p.fInv, p.zzp_context, p.factors, p.e);
    p.sigmaX = getsigmaXVal(p.wLen, p.mBar, p.maxFactor, p.r);

    seekPastChar(s, ']');  // this function is defined in tools.cpp
    return s;
}

// binary I/O - write all the variables to the file. The handle of the open file is an input to the function The function returns the number of items written
long TDMatrixParams::writeToFile(FILE* handle) const
{
    FHE_TIMER_START;
    #ifdef DEBUGPRINT
    char cwd[1024];
    getcwd(cwd, sizeof(cwd));
    cout << "write TDMatrixParams, current directory = " << cwd << endl;
    #endif

    long count = fwrite(&n, sizeof(n), 1, handle);
    count += fwrite(&m, sizeof(m), 1, handle);
    count += fwrite(&mBar, sizeof(mBar), 1, handle);
    count += fwrite(&e, sizeof(e), 1, handle);
    count += fwrite(&r, sizeof(r), 1, handle);
    count += fwrite(&kFactors, sizeof(kFactors), 1, handle);
    count += fwrite(factors.elts(), sizeof(long), factors.length(), handle);

    return count;
}

//Binary IO - reads the class parameters from the file. The open file handle
//is provided as input to the function, and the number of items read is
//returned by the function.
long TDMatrixParams::readFromFile(FILE* handle)
{
    FHE_TIMER_START;
    long count = 0;

#ifdef DEBUGPRINT
    char cwd[1024];
    getcwd(cwd, sizeof(cwd));
    cout << "read TDMatrixParams from current directory = " << cwd << endl;
#endif

    count += fread(&n, sizeof(n), 1, handle);
    count += fread(&m, sizeof(m), 1, handle);
    count += fread(&mBar, sizeof(mBar), 1, handle);
    count += fread(&e, sizeof(e), 1, handle);
    count += fread(&r, sizeof(r), 1, handle);
    count += fread(&kFactors, sizeof(kFactors), 1, handle);

    factors.SetLength(kFactors);
    count += fread(factors.elts(), sizeof(long), factors.length(), handle);

    //initialize stash
    stash.SetLength(kFactors);
    for (long i = 0; i < kFactors; i++)
      stash[i].init(factors[i]);

    wLen = n * kFactors * e;

    maxFactor = setFactors(fInv, zzp_context, factors, e);
    sigmaX = getsigmaXVal(wLen, mBar, maxFactor, r);

    return count;
}

//check if all the variables in the two classes are equal if Yes, return true. Else, false.

bool operator==(const TDMatrixParams& p, const TDMatrixParams& q)
{
    return (p.n == q.n && p.m == q.m && p.mBar == q.mBar && p.e == q.e
            && p.r == q.r && p.wLen == q.wLen && p.kFactors == q.kFactors
            && p.factors == q.factors);
}

#if 0
// Deprecated, kept for debuging purposes
TDMatrixParams::TDMatrixParams(vec_l &vFactors, long lm,long rSigma,long llogQ,long lq, long nIn, long eL)
{
    FHE_TIMER_START;
    factors = vFactors;
    r = rSigma;
    kFactors = llogQ;
    n = nIn;
    m = lm;//SMS.NumCols();
    mBar = m - n*(kFactors*eL);
    wLen = n*kFactors*eL;
    e = eL;

    // For each factor f_i, initialize zz_p context for F_i=f_i^e,
    // and also compute f_i^{-1} mod F_j for all i<j
    maxFactor = setFactors(fInv, zzp_context, factors, e);

    sigmaX = getsigmaXVal(wLen, mBar, maxFactor, r);
    // The Gaussian parameter from which we can sample with a trapdoor

    stash.SetLength(kFactors);
    for (long i = 0; i < kFactors; i++)
      stash[i].init(factors[i]);
}
#endif
