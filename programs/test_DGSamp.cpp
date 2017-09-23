/* Copyright (C) 2017 IBM Corp.
 * Licensed under the Apache License, Version 2.0 (the "License"); 
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
test_DGSamp: a program for testing the DGSamp module
********************************************************************/

// ===> This test program is obsolete, need to be updated <===

// test_DGSamp.cpp - Test trapdoor sampling, including statistical
// tests to ensure that we get the right probability distribution.

#include <cassert>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stack>
#include <string>
#include <vector>

// #include <mpfr.h>

#include <NTL/mat_ZZ_p.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/matrix.h>
#include <NTL/vec_vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include "../mat_l.h"
NTL_CLIENT


#include "../vec_l.h"
#include "../DGaussSampler.h"
#include "../utils/timing.h"
#include "../TDMatrix.h"
#include "../utils/argmap.h"

#include <stdio.h>
#include <string.h>

bool MovetoNextxTest(vec_l &xTest, int VecLen, int noVars);

int TestProbabilities(TDMatrixParams& TDMP, long nVars);

int main(int argc, char *argv[])
{
    ArgMapping amap;

    long k = 2;
    amap.arg("k", k, "number of factors (~ log q)");

    long n = 2;
    amap.arg("n", n, "the small dimension"); // no default info

    long e=2;
    amap.arg("e", e, "the power of each factor");

    long m=0;
    amap.arg("m", m, "the dimension m", "0: set as n*(2+k*e)");

    long nVars=2;
    amap.arg("nVars", nVars, "# of values to check on each side of the mean (in addition to the mean itself), Total vals checked = 2*nVars+1");

    amap.parse(argc, argv); // parses and overrides initail values
    zz_pPush push; // backup NTL's current modulus
    TDMatrixParams TDMP(n,k,e,m);

    TestProbabilities(TDMP, nVars);
    cout << "maxSigma = "<< Gaussian1Dsampler::maxSigma << endl;
    cout << "maxSample = "<< Gaussian1Dsampler::maxSample << endl << endl;

    printAllTimers();

#ifdef CodeBlocks
    cin.get();
#endif

    return 0;
}

/*
Test that the probabilities of getting each sample are correct
*/

int TestProbabilities(TDMatrixParams& TDMP, long nVars)
{
#ifdef DEBUG
    NTL::SetSeed((unsigned char*)"TestProb",8);
#endif
    TDMP.m = 1; //1 dimensional
    //Test the probability of getting each sample
    DiscreteGaussianSampler DGS;
    //test the DGaussSampler

    long VectorLength = TDMP.m;

    Vec<double> zeroVector(INIT_SIZE, TDMP.m);
    clear(zeroVector);

    vec_l Counter;
    int noVars = 2*nVars+1; //5
    long counterSize = pow(noVars,VectorLength);
    int HalfVars = (noVars-1)/2;
    int HalfCounterSize = (counterSize - 1)/2;
    Counter.SetLength(counterSize);//NoVars);

    for (long i = 0; i < counterSize; i++)
    {
        Counter.put(i,0);
    }
#ifdef DEBUG
    cout << "Counter=" << Counter << "\n";
#endif

    mat_l SMS;
    SMS.SetDims(TDMP.m,TDMP.m); //SMS larger or equal to RI *(2+sigma G)*RTI

    mat_l A;
    A.SetDims(TDMP.m,TDMP.m);
    ident(A,TDMP.m);

    //A.put(1,0,1);

    mat_l ATrans;
    square(SMS,A);
    //mul(SMS,A,ATrans);
    mat_l SMSinv;

    if (TDMP.m==1)
    {
        SMS*=16;
        ident(SMSinv,1);
    }
    else
    {
        ZZ zz_det;
        mat_ZZ zz_SMS;
        zz_SMS = conv<mat_ZZ>(SMS);
        mat_ZZ zz_SMSINV;
        inv(zz_SMSINV,zz_SMS);
        determinant(zz_det,zz_SMS);

#ifdef DEBUG
        cout << "SMS= " << zz_SMS << "\n";
        cout << "zz_SMSINV= " << zz_SMSINV << "\n";
#endif

        SMSinv=conv<mat_l>(zz_SMSINV);
    }

#ifdef DEBUG
    {
        cout << "SMS = \n";
        cout << SMS << "\n";
        cout << SMSinv << det << "\n";
    }
#endif

    vec_l xOut;

    zeroVector.SetLength(TDMP.m);
    for (long i = 0; i< TDMP.m; i++)
        zeroVector.put(i,0);

    DGS.InitSampler(SMS);

    long NoTrials = 0;
    long minMidVal = 1000000;//10000*(nVars^2);
    //long var = zeroVector.get(0);

    while (Counter[HalfCounterSize] < minMidVal)
    {
        xOut.SetLength(0);

        DGS.SampleDiscreteGaussian(xOut,zeroVector);
        NoTrials++;

        long iIndex = 0;
        for (long j=0; j< TDMP.m; j++)
            iIndex += xOut[j] * (pow(noVars,j));

        //long iIndex = xOut[0];
        if (abs(iIndex) > HalfCounterSize) //index out of range
            continue;

        Counter[iIndex + HalfCounterSize]+=1;

#ifdef DEBUG
        long var2 = Counter.get(0);//iIndex+HalfVars);d
        cout << var2;
        var2 = Counter.get(1);
        cout << ", " << var2;
        var2 = Counter.get(2);
        cout << ", " << var2 << "\n";
#endif // DEBUG

        // if (Counter[HalfCounterSize] > minMidVal)
        //     bReachedLimit = true;

    }

#ifdef DEBUG
    cout << "\n";

    cout << "Counter=" << "\n";

    for (long ii = ((counterSize-1)/2 -HalfVars); ii <= ((counterSize-1)/2 +HalfVars); ii++)
        cout << "Counter["<<ii<<"]=" << Counter[ii] << "\n";

    cout << "NoTrials = " << NoTrials << "\n";
#endif // DEBUG

    //calculate the error value
    int counterIndex;
    double probVal;
    Vec<double> Prob_vec;
    vec_l val_vec;
    Vec<double> error_vec;

    //calculate avg error values for all results
    //prob = a * exp(-0.5 * X * SigInv * Xt);
    vec_l xTest;
    vec_l xTestNorm;
    xTest.SetLength(TDMP.m);
    xTestNorm.SetLength(TDMP.m);
    Prob_vec.SetLength(counterSize);
    val_vec.SetLength(counterSize);

    for (long ii = 0; ii < TDMP.m; ii++)
        xTest[ii] = 0;

    error_vec.SetLength(counterSize);
    //test for each combination

    while (1)
    {
        //generate counterindex
        counterIndex = 0;
        for (long i = 0; i < VectorLength; i++)
        {
            counterIndex += xTest[i];
            if (i < VectorLength-1)
                counterIndex*=noVars;
        }

        //p = a*exp(-0.5* X * SigInf * XT);

        xTestNorm = xTest - HalfVars;

        double sigma = sqrt(SMS.get(0,0));

        //calculate the probability value of each variable

        double lowerHalf = (double) xTestNorm[0] - 0.5;
        double upperHalf = (double) xTestNorm[0] + 0.5;

#ifdef DEBUG
        cout << "lowerHalf = " << lowerHalf << ", upperHalf = " << upperHalf  << "\n";
#endif

        //tested mean = 0
        lowerHalf/=sigma;
        lowerHalf/=sqrt(2);
        upperHalf/=sigma;
        upperHalf/=sqrt(2);
        double erfUpper = erf(upperHalf);
        double erfLower = erf(lowerHalf);
        probVal = erfUpper - erfLower;
        probVal*=0.5;
#ifdef DEBUG
        cout << "lowerHalf = " << lowerHalf << ", upperHalf = " << upperHalf << ", probVal = " << probVal << "\n";
#endif

        Prob_vec[counterIndex] = probVal;
        probVal*= ((double)NoTrials/Counter.get(counterIndex));
        error_vec.put(counterIndex,probVal);

        if (!MovetoNextxTest(xTest, TDMP.m, noVars))
            break;
    }

#ifdef DEBUG

    cout << "Prob_mat = \n";
    for (counterIndex = 0; counterIndex < counterSize; counterIndex++)
    {
        double dVal = Prob_vec[counterIndex];

        if (dVal != -1)
        {
            cout << "Index = " << counterIndex << "\n";
            cout << dVal << " ";
            dVal = error_vec[counterIndex];
            cout << " error_vec =" << dVal << "\n";
        }


    }
    cout << "\n";
#endif


    for (counterIndex = 0; counterIndex < counterSize; counterIndex++)
    {
        double dVal = Prob_vec[counterIndex];

        if (dVal != -1)
        {
            dVal = error_vec[counterIndex];

            if ((dVal > 0.98) && (dVal < 1.02))
            {
                cout << "PASSED" << counterIndex << ", difference =  " << (abs((dVal -1) * 100)) << "%" << endl;
            }
            else
            {
                cout << "FAILED " << counterIndex << ", difference = " << (abs((dVal -1) * 100)) << "%" << endl;

            }
        }
    }

    return 0;
}

//increase xTest vector in lexicographical order
//returns true as long as moves to next vector, false when it is done

bool MovetoNextxTest(vec_l &xTest, int VecLen, int noVars)
{
    if (xTest[VecLen-1]==(noVars-1))
        return false; //done

    for (long i = 0; i < VecLen; i++)
    {
        xTest[i] += 1;
        if (xTest[i] >= noVars)
            xTest[i] = 0;
        else
            break;
    }
    return true;

}
