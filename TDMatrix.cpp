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
 * TDMatrix.cpp - Matrixes with trapdoors,
 *   implementing the Efficient algorithm SampleD - Algorithm 3 from MP12.
 *****************************************************************************/
#include <stdexcept>
#include <ctime>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_lzz_p.h>
#include <iostream>
#include <fstream>

NTL_CLIENT
#include "utils/timing.h"
#include "TDMatrix.h"
#include "mat_l.h"

//#define DEBUGPRINT
//#define DEBUG
//#define PRINTDOT //show how program progresses

// A static variable for debugging purposes
std::atomic<int> TDMatrix::maxSample(0); // keep the largest sample ever drawn

/***** Some local functions *****/

// Compute the covariance matrix sigmaP = sigmaX*I - sigmaG*[R/I]*(R^t|I)
static void
computeCovarianceMatrix(mat_l& covMat,const mat_l& R,long sigmaG,long sigmaX);

// Multiply a vector by the G matrix modulo the current NTL modulus, y=G*x.
// G is and n-by-m matrix with m=n*e*numOfFactors, and is defind by means
// of the vector of factors and the exponent e, each row of G is mostly zero,
// except for a progression
//    ( 1, f1,...,f1^e,  (f1^e f2), (f1^e f2^2), ..., (f1^e...fk^{e-1}) )
template<class T>
void multByG(vec_zz_p &y, const Vec<T> &x,
             long n, const vec_l& factors, long e)
{
    FHE_TIMER_START;
    long k = factors.length();
    long m = n*k*e;

    assert(x.length()==m);
    y.SetLength(n);

    for (long row=0; row<n; row++)   // Multiply x by next row of G
    {
        long index = row*k*e;      // first non-zero entry in this row
        zz_p val = to_zz_p(1L);  // value in current entry of G
        y[row] = to_zz_p(0L);    // accumulated sum
        for (long f=0; f<k; f++)     // go over all the factors
        {
            zz_p fact = to_zz_p(factors[f]); // convert to zz_p
            for (long i=0; i<e; i++)    // use each factor e times
            {
                zz_p addedTerm = val;
                addedTerm *= x[index++];  // update sum, advance index
                y[row] += addedTerm;
                val *= fact;              // compute next entry in this row of G
            }
        } // end of current factor
    } // end of current row
}
// Explicit instantiations for Vec<long> and Vec<zz_p>
template void multByG<zz_p>(vec_zz_p &y, const Vec<zz_p> &x,
                            long n, const vec_l& factors, long e);
template void multByG<long>(vec_zz_p &y, const Vec<long> &x,
                            long n, const vec_l& factors, long e);


// Choose a (pseudo)random A with a trapdoor R s.t. A*[R/I]=G
// A=(aBar| aPrime ) with aPrime=G-aBar*R isn't computed explicitly here
void TDMatrix::initTDmatrix(TDMatrixParams& prms, const CRTmatrix* abar)
{
    FHE_TIMER_START;
    params = &prms;

    if (abar)   // aBar is specified by the caller
        aBarCRT = *abar;
    else        // Choose aBar at random
        randomFill(aBarCRT, prms.n, prms.mBar, prms);

    mat_l sigmaP;
    do
    {
        // Choose R as a random small matrix
        setSmall(trapDoorMatR, params->mBar,params->wLen, params->r);

        // Compute the covariance matrix sigmaP = sigmaX*I - sigmaG*[R/I]*(R^t|I)

        computeCovarianceMatrix(sigmaP, trapDoorMatR,
                                /*sigmaG=*/prms.r *prms.maxFactor, prms.sigmaX);

#ifdef DEBUGPRINT
        long minVal =0;
        long maxVal=0;
        for (long i = 0; i < sigmaP.NumCols(); i++)
            for (long j = 0; j < sigmaP.NumRows(); j++)
            {
                if (sigmaP[i][j] < minVal)
                    minVal = sigmaP[i][j];
                if (sigmaP[i][j] > maxVal)
                    maxVal = sigmaP[i][j];
             }

        ZZ q = prms.getQ();
        cout << "minVal=" << minVal << endl;
        cout << "maxVal=" << maxVal << endl;
        cout << "q=" << q << endl;
        cout << "after computeCovarianceMatrix"  << endl;
#endif // DEBUGPRINT
    }
    while (!(gaussSamp.InitSampler(sigmaP)));  // error, try again

#ifdef DEBUGPRINT
        cout << "after gaussSamp.InitSampler"  << endl;
#endif // DEBUGPRINT
}


// Compute the covariance matrix sigmaP = sigmaX*I - sigmaG*[R/I]*(R^t|I)
// Note that [R/I]*(R^t|I) is a block matrix with this form:
//    RR = ( R^2 | R^t )
//         (  R  |  I  )
void computeCovarianceMatrix(mat_l& covMat, const mat_l& R, long sigmaG,long sigmaX)
{

    FHE_TIMER_START;

    long mBar = R.NumRows();
    long wLen = R.NumCols();
    long m = mBar + wLen;

    covMat.SetDims(m,m); // allocate space

    mat_l tmp(INIT_SIZE, mBar, mBar);
    square(tmp, R);   // Set the top-left mBar-by-mBar as R*R^t

    // Reset the top-left to sigmaX*I - sigmaG*R*R^t
    for (long i=0; i<mBar; i++) for (long j=0; j<mBar; j++)
        {
            if (i==j)
                covMat[i][j] = sigmaX - sigmaG*tmp[i][j];
            else
                covMat[i][j] = -sigmaG*tmp[i][j];
        }

    // Set the top-right ro -sigmaG*R^t and the bootom-left to -sigmaG*R
    for (long i=0; i<mBar; i++) for (long j=0; j<wLen; j++)
        {
            covMat[i][j+mBar] = covMat[j+mBar][i] = -sigmaG * R[i][j];
        }

    // Finally set the bottom-right to (sigmaX-sigmaG)*I;
    for (long i=mBar; i<m; i++) for (long j=mBar; j<m; j++)
        {
            if (i==j) covMat[i][j] = sigmaX -sigmaG;
            else      covMat[i][j] = 0;
        }
}


// Compute A explicitly, A = ( aBar | G - aBar*R )
void TDMatrix::getA(CRTmatrix& A) const
{
    FHE_TIMER_START;
    A.params = params;

    long n = params->n;
    long m = params->m;
    long kFactors = params->kFactors;
    long mBar = params->mBar;

    // Compute A modulo each factor separately

    zz_pPush ppush; // backup NTL's current modulus
    A.SetLength(kFactors);
    for (long f=0; f<kFactors; f++)
    {
        params->zzp_context[f].restore();
        mat_zz_p& aMod = A[f];
        const mat_zz_p& aBarMod = aBarCRT[f];

        aMod.SetDims(n,m); // allocate space

        // Copy aBar
        for (long i=0; i<n; i++) for (long j=0; j<mBar; j++)
            {
                aMod[i][j] = aBarMod[i][j];
            }

        // Compute aBar * R
        mat_zz_p tmp;
        mat_zz_p rMod;
        conv(rMod, trapDoorMatR); // conver to zz_p format
        mul(tmp, aBarMod, rMod); // tmp = aBar * R (mod P_f)

        // Copy -aBar*R to left part of A
        for (long i=0; i<n; i++) for (long j=mBar; j<m; j++)
            {
                aMod[i][j] = -tmp[i][j-mBar];
            }

        // Add G to left part of A
        for (long i=0; i<n; i++)    // one row at a time
        {
            // Add (1, f1, f1^2, ..., f1^e...fk^{e-1}) to each row

            long index = (i * kFactors * params->e) + mBar; // 1st index to add to

            zz_p val = to_zz_p(1L);      // The value to add to the next entry
            for (long j=0; j < kFactors; j++)
                for (long ei=0; ei<params->e; ei++, index++)
                {
                    aMod[i][index] += val;
                    val *= params->factors[j];// Change value to add before next entry
                }
        }
    }

}

/** sampleWithTrapdoor - Implementation of Efficient algorithm SampleD -
 * Algorithm 3 from MP12. Choose p at random and z to match the needed
 * solution. This debugging function also returns p and z.
 **/
int TDMatrix::sampleWithTrapdoor(vec_l &xOut, vec_l& p, vec_l& z,
                                 Vec<vec_zz_p>& U) const
{
    FHE_TIMER_START;
    long mBar = params->mBar;
    long m = params->m;
    long n = params->n;
    long e = params->e;
    long kFactors = params->kFactors;

    //choose perturbation P from (0,sigma)

    Vec<double> zeroVector(INIT_SIZE, m);
    clear(zeroVector);

    gaussSamp.SampleDiscreteGaussian(p, zeroVector);
#ifdef DEBUGPRINT
    cout << "\np = " << p << endl;
#endif
    // Split the pertubation vector in two
    vec_l p1(INIT_SIZE, mBar);
    for (long i = 0; i < mBar; i++)
        p1[i] = p[i];

    vec_l p2(INIT_SIZE, params->wLen);
    for (long i = mBar; i < m; i++)
        p2[i - mBar] = p[i];

    //compute wBar and w

    /* We need to compute v = u-A*p = u-wBar-w (in CRT representation), where
     *
     *    A*p = (Abar| Aprime ) * [p1/p2] = Abar*p1 + Aprime*p2
     *        = Abar*p1 + (G-Abar*R)*p2   = Abar*(p1-R*p2) +  G*p2
     *                                      \___wBar____/    \_w_/
     */

    // p1-R*p2 is computed over the integers, stored in p1
    vec_l tmp(INIT_SIZE, mBar);
    mul(tmp, trapDoorMatR, p2); // R*p2
    p1 -= tmp;                  // p1-R*p2

    // All other vectors are computed modulo each of the factors

    // FIXME: There's room for parallelization here, but in our case we have
    // thousands of calls to sampleWithTrapdoor that can be run in parallel,
    // so there is no reason to parallelize also at this lower level.

    zz_pPush push; // backup NTL's current modulus

    Vec<vec_zz_p> vVec = U; // allocate space and initialize v=u
    vec_zz_p wMod, wBarMod;
    for (long iE = 0; iE < kFactors; iE++)
    {
        params->zzp_context[iE].restore();

        // Compute Abar*(p1-R*p2) modulo the current factor
        const mat_zz_p& aBarMod = aBarCRT[iE];
        mul(wBarMod, aBarMod, conv<vec_zz_p>(p1));  // wBar = Abar*(p1-R*p2)
        // convert p1-R*p2 to vec_zz_p and multiply by aBar modulo each factor

        vVec[iE] -= wBarMod; // u - wBar

        multByG(wMod, p2, n, params->factors, e); // w = G*p2
        vVec[iE] -= wMod;    // u - wBar - w
    }

    // Now that we computed v modulo all the primes, sample z s.t. G*z=v

    sampleG(z, vVec, params);
#ifdef DEBUGPRINT
    cout << "z = " << z << endl;
#endif
    // finally, return x = p + [R/I]*z = p + [R*z / z]

    xOut = p;
    mul(tmp, trapDoorMatR, z); // add R*z to top mBar entries
#ifdef DEBUGPRINT
    cout << "R*z = " << tmp << endl;
#endif
    for (long i=0; i<mBar; i++)
    {
      xOut[i] += tmp[i];
      long absX = abs(xOut[i]);
      if (absX > TDMatrix::maxSample) TDMatrix::maxSample.store(absX);
         // not quite thread-safe, the stored value may not be the maximum
    }

    // add z to bottom m-mBar entries
    for (long i = mBar; i< params->m; i++)
    {
      xOut[i] += z[i - mBar];
      long absX = abs(xOut[i]);
      if (absX > TDMatrix::maxSample) TDMatrix::maxSample.store(absX);
         // not quite thread-safe, the stored value may not be the maximum
    }
    return 0;
}

//Write TDMatrix different parameters to file. The handle to the open
//file is sent as input to the function. returns number of items read by
//the function
long TDMatrix::writeToFile(FILE* handle)
{
    FHE_TIMER_START;
    long count = 0;

    count += params->writeToFile(handle); // write the params
    count += aBarCRT.writeToFile(handle); // write the Abar matrix

    long numRows = trapDoorMatR.NumRows();
    long numCols = trapDoorMatR.NumCols();
    count += fwrite(&numRows,sizeof(long),1,handle);
    count += fwrite(&numCols,sizeof(long),1,handle);
    for (long i = 0; i < numRows; i++)  // write the rows of R, one at a time
    {
        count+= fwrite(trapDoorMatR[i].elts(), sizeof(long), numCols, handle);
        // each row is implemented as a C vector
    }

    count += gaussSamp.writeToFile(handle); // write the DGaussSampler

    return count;
}

//Read the TDMatrix parameters from a flie. The handle to the open file is
//received as input to this function, as well as a pointer to the TDMatrixParams
//if the pointer is not provided and params as not been initialized before, an
//error occurs, as a pointer to this structure is needed
//returns number of items read by the function

long TDMatrix::readFromFile(FILE* handle, TDMatrixParams* prmBuf)
{
    FHE_TIMER_START;
    assert(params != NULL || prmBuf != NULL); // some pointer must be provided
    long count=0;

    // If buffer is given, make params point to it and don't overwrite from iput
    if (prmBuf != NULL)
    {
        TDMatrixParams p;
        count += p.readFromFile(handle);
        assert(p == *prmBuf); // sanity check
        params = prmBuf;      // point to given params
    }
    else
        count += params->readFromFile(handle); // overwrite params from input

    count += aBarCRT.readFromFile(handle,prmBuf); // read the Abar matrix

    long numRows,numCols;
    count += fread(&numRows,sizeof(long),1,handle);
    count += fread(&numCols,sizeof(long),1,handle);
    trapDoorMatR.SetDims(numRows,numCols);
    for (long i = 0; i < numRows; i++)  // write the rows of R, one at a time
    {
        count+= fread(trapDoorMatR[i].elts(), sizeof(long), numCols, handle);
        // each row is implemented as a C vector
    }

    count += gaussSamp.readFromFile(handle); // read the DGaussSampler

    return count;
}

// operator==, //compares all variables in TDMatrix
bool operator==(const TDMatrix& A, const TDMatrix& B)
{
    FHE_TIMER_START;
    if (A.getParams() != B.getParams())
        return false;

    if (A.getABar() != B.getABar())
        return false;

    if (A.getR() != B.getR())
        return false;

    return true;
}

#ifdef DEBUG
static void printTrapdoor(const mat_l& td, const TDMatrixParams& params)
{
    FHE_TIMER_START;
    double maxR = 0;
    double AvgR=1;
    //  long TotalItems;

    for (long i=0; i < params.mBar; i++)
        for (long j=0; j < params.wLen; j++)
        {
            if (abs(td[i][j]) > maxR)
                maxR = abs(td[i][j]);
            AvgR *= (double)(i*params.wLen+j)/(i*params.wLen+j+1);
            AvgR+=(double)abs(td[i][j]) /(i*params.wLen+j+1);
        }
    cout << "maxR = " << maxR << ", AvgR = " << AvgR << "\n";

    mat_l rSquare;
    //how large is R^2?
    square(rSquare,td);
    //  cout << "trapDoorMatR=" << td << "rSquare= " << rSquare << "\n";

    maxR = 0;
    AvgR=1;
    for (long i=0; i < params.mBar; i++)
        for (long j=0; j < params.mBar; j++)
        {
            if (abs(rSquare[i][j]) > maxR)
                maxR = abs(rSquare[i][j]);
            AvgR *= (double)(i*params.wLen+j)/(i*params.wLen+j+1);
            AvgR+=(double)abs(rSquare[i][j]) /(i*params.wLen+j+1);
        }
    cout << "maxRSquare = " << maxR << ", AvgRSquare = " << AvgR << "\n";
}
#endif //DEBUG
