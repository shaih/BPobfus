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
/**********************************************************
 * test_TDMatrix: a program for testing the TDMatrix module
 * prints out PASS/FAIL information. If a certain bit fails,
 * the bit index is returned in the message
***********************************************************/
//#include <mpfr.h>
#include <NTL/lzz_p.h>

NTL_CLIENT
#include "../utils/argmap.h"
#include "../TDMatrix.h"

//#define VERBOSE //check extra information?

int main(int argc, char *argv[])
{
#ifdef DEBUG
    NTL::SetSeed((unsigned char*)"test_sampleG",12);
#endif
    ArgMapping amap;
    long k = 2;
    amap.arg("k", k, "number of factors (~ log q)");
    long n = 2;
    amap.arg("n", n, "the small dimension"); // no default info
    long e=2;
    amap.arg("e", e, "the power of each factor");
    long m=0;
    amap.arg("m", m, "the dimension m", "0: set as n*(2+k*e)");

    amap.parse(argc, argv); // parses and overrides initail values
    // for each parameter, one per line,

    if (m < n*(2+k*e)) m = n*(2+k*e);

    TDMatrixParams params(n,k,e,m);
    //    Vec<long>& factors = params.factors;

    cout << "n=" << n << ", m=" << m << ", k="
         << k     << ", e=" << e
	 << ", sigmaX=" << params.sigmaX << ", q=" << params.getQ() << endl;

    long mBar = params.mBar;
    TDMatrix aMat(params);
    CRTmatrix A;
    aMat.getA(A);

    
    const mat_l& R = aMat.getR();
    #ifdef VERBOSE
    if (R.NumRows()<50 && R.NumCols()<50)
        cout << "\n R=" << R << endl;
    #endif

    mat_l RI(INIT_SIZE, m, n*k*e);
    for (long i=0; i<mBar; i++) for (long j=0; j<n*k*e; j++)
        {
            RI[i][j] = R[i][j];
        }
    for (long i=0; i<n*k*e; i++) for (long j=0; j<n*k*e; j++)
        {
            if (i==j) RI[i+mBar][j] = 1;
            else      RI[i+mBar][j] = 0;
        }

    #ifdef VERBOSE
    for (long i=0; i<k; i++)
    {
        params.zzp_context[i].restore();
        if (params.m < 50)
            cout << "\n A mod " << zz_p::modulus() << "=" << A[i] << endl;

        mat_zz_p A_RI;
        mul(A_RI, A[i], conv<mat_zz_p>(RI));
        if (params.m < 50)
            cout << "\n A*[R/I] mod " << zz_p::modulus() << "="
                 << A_RI << endl;
    }
    #endif

    // Compute a random syndrome
    Vec<vec_zz_p> syndrome(INIT_SIZE, k);
    for (long i=0; i<k; i++)
    {
        params.zzp_context[i].restore();
        syndrome[i].SetLength(n);
        for (long j=0; j<n; j++) syndrome[i][j] = random_zz_p();

        #ifdef VERBOSE
        if (params.m < 50)
            cout << "syndrome mod "
                 << power_long(factors[i],e)<<" ((factor "<<i<<")^"<<e+1<<")= "<<syndrome[i]<<endl <<flush;
        #endif
    }

    // Sample x such that A*x = syndrome
    Vec<long> x;
    aMat.sampleWithTrapdoor(x, syndrome);
    // vec_l p, z;
    // aMat.sampleWithTrapdoor(x, p, z, syndrome); // debugging version

    #ifdef VERBOSE
    if (params.m < 50) cout << "sampled vector="<<x<<endl;
    #endif

    assert(x.length()== m);
    // check that we have the right answer modulo p^e for all factors

    bool bSuccess = true;

    for (long i=0; i<k; i++)
    {
        params.zzp_context[i].restore();
        vec_zz_p xMod = conv<vec_zz_p>(x);
        vec_zz_p uu;
        mul(uu, A[i], xMod);
        if (syndrome[i]!=uu)
        {
            cout << "bit " << i << " failed!" << endl;
            bSuccess = false;

            //assert(syndrome[i]==uu);
        }
    }
    if (bSuccess == true)
        cout << "\nPASSED\n";
    cout << "maxSigma="<<Gaussian1Dsampler::maxSigma
         << ", maxSample="<<TDMatrix::maxSample << endl;

#ifdef DEBUG
    printAllTimers(cout);
#endif

#ifdef CodeBlocks
    cin.get();
#endif

    return 0;


}
