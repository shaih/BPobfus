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
/*****************************************************************
 * test_CRTmatrix: a program for testing the CRTmatrix module
 *****************************************************************/
//#include <mpfr.h>
#include <NTL/lzz_p.h>

NTL_CLIENT
#include "../utils/argmap.h"
#include "../TDMatrix.h"

#if defined(__unix__) || defined(__unix) || defined(unix)
#include <sys/time.h>
#include <sys/resource.h>
#endif

#include <NTL/BasicThreadPool.h>

#define NFACTORS 8

long tinyFactors[NFACTORS] = { 2, 3, 5, 7, 11, 13, 17, 19 };

int main(int argc, char *argv[])
{
#ifdef DEBUG
    NTL::SetSeed((unsigned char*)"test_CRTmatrix",15);
#endif
    ArgMapping amap;
    long k = 3; //2
    amap.arg("k", k, "number of factors (~ log q)");
    long n = 5; //2
    amap.arg("n", n, "the small dimension");
    long e=3; //2
    amap.arg("e", e, "the power of each factor");
    long m=100; //0
    amap.arg("m", m, "the dimension m", "0: set as n*(2+k*e)");
    long threads = 2;
    amap.arg("threads", threads, "how many threads to run");

    amap.parse(argc, argv); // parses and overrides initail values
    // for each parameter, one per line,

    TDMatrixParams prms(n,k,e,m);
    m = prms.m;
    cout << "n="<<n<<", m="<<m<<", k="<<k<<", e="<<e<<", q="<<prms.getQ()
         << endl << endl;

    // Setup a htread pool with almost as many threads as there are cores
    if (threads<1) threads=1;
    else if (threads > numCores()-2) { // numCores() defined in utils/tools.h
      threads = numCores()-2;
      if (threads<1) threads=1;
    }
    NTL::SetNumThreads(threads);
    cout << "Using "<<threads<< " threads\n";

    CRTmatrix A(&prms);
    CRTmatrix B(&prms);
    CRTmatrix C(&prms);
    CRTmatrix I(&prms);
    CRTmatrix T1(&prms);
    CRTmatrix T2(&prms);
    CRTmatrix T3(&prms);
    CRTmatrix T4(&prms);

    generateMultiPair(A, B, prms.m);

    generateMatrixPair(A, B, n);
    C.invert(B);
    bool bTestRes = A==C;
    if (bTestRes==false)
    {
        cout << "Failed Generate Matrix Pair" << endl;
    }
    else    
        cout << "Passed Generate Matrix Pair" << endl;

    bTestRes = A==A;

    if (bTestRes==false)
    {
        cout << "Failed Equality (==) Test" << endl;
    }
    else
         cout << "Passed Equality (==) Test" << endl;

    I.identity(n);

    //identity test

    bool bFailedIDTest = false;
    bool bFailedColumnTest = false;

    for (long lFactor=0; lFactor < prms.kFactors; lFactor++)
    {
        const mat_zz_p& aMod = I[lFactor];
        prms.zzp_context[lFactor].restore();

        for (long i = 0; i < n; i++)
        {
            CRTvector crtCol;
            I.getColumn(crtCol, i);
            for (long j = 0; j < n; j++)
            {
                if (i==j)
                {
                    if (aMod[i][j]!=1)
                    {
                        cout << "Failed Identity Test" << endl;
                        bFailedIDTest = true;
                    }
                    if (!IsOne(crtCol[lFactor][j]))
                    {
                        cout << "Failed Get Column Test" << endl;
                        bFailedColumnTest = true;
                    }
                }
                else
                {
                    if (aMod[i][j]!=0)
                    {
                        cout << "Failed Identity Test" << endl;
                        bFailedIDTest = true;
                    }
                    if (!IsZero(crtCol[lFactor][j]))
                    {
                        cout << "Failed Get Column Test" << endl;
                        bFailedColumnTest = true;
                    }
                }
            }
        }
    }

    if (bFailedIDTest == false)
        cout << "Passed Identity Test" << endl;
    if (bFailedColumnTest == false)
        cout << "Passed Get Column Test" << endl;

//cout << "A="<<A<<endl;
//cout << "I =" << I << endl;

// sanity-check
//    assert(A != I);

    A *= B;
    if (A == I) cout << "PASSED generate Matrix Pair\n";
    else        cout << "FAILED generate Matrix Pair\n";

    if (!(T1.NumCols() ==n))   //test NumCols()
        bTestRes = false;

    if (!(T1.NumRows() ==n))   //test NumRowls
        bTestRes = false;

    if (bTestRes==false)
        cout << "FAILED rows/cols num checking" << endl;
    else cout << "PASSED rows/cols num checking" << endl;


    mat_l longA;
    do
    {
        setSmall(longA, prms.n, prms.n, prms.r);
    }
    while (!B.invert(longA)); // Repeat until you get an invertible A mod q
    //cout << "longA=" << longA << "\n" << "B=" <<B << endl;

    A = B;
    
    B.leftMultBy(longA); //test left multiply by
    A *= longA; //test the *= multiplication


    if (A==I && B==I) cout << "PASSED create invertible matrix\n";
    else              cout << "FAILED create invertible matrix\n";

    T1.randomFill(prms.n,prms.n);
    //cout << "T1=" << T1 << endl;
    T2.randomFill(prms.n,prms.n);
    //cout << "T2=" << T2 << endl;
    // bool bIsSmall= T1.isSmall(prms.r);
    // bIsSmall= B.isSmall(prms.r);

    Vec<vec_zz_p> TestCol;
    B.getColumn(TestCol,1);
    //cout << TestCol << "\n";

    T3= T1;
    //cout << "T3=" << T3 << endl;
    //cout << "T1=" << T1 << endl;
    T1-=I;          //test -= multiplication
    T3-=I;

    //cout << "I=" << I << endl;
    //cout << "T1=" << T1 << endl;

    bool bPassedPMtests = true;

    if (!(T3==T1))
    {
        cout << "FAILED -= test" << endl;
        bPassedPMtests = false;
    }
    T1+=I;          //test += multiplication
    T3+=I;

    if (T3!=T1)
    {
        cout << "FAILED += test" << endl;
        bPassedPMtests = false;
    }
    //cout << "T1=" << T1 << ", T2=" << T2 << "\n";
    T1-=T2;
    //cout << "T1=" << T1 << "\n";
    T1+=T2;

    if (T3!=T1)
    {
        cout << "FAILED +=, -= tests" << endl;
        bPassedPMtests = false;
    }
    if (bPassedPMtests)
         cout << "PASSED +=, -= tests" << endl;


    FILE* handle = fopen("CRMatrixTest.dat", "wb");
    T2.writeToFile(handle);
    fclose(handle);

    handle = fopen("CRMatrixTest.dat", "rb");
    T4.readFromFile(handle,&prms);

    if (bTestRes)
        cout << "PASSED function test" << endl;
    else
        cout << "FAILED functions test" << endl;

   // cout << "T2=" << T2 << ", T4=" << T4 << endl;
    if (T2==T4)
        cout << "Read/Write Passed!" << endl;
    else
        cout << "Read/Write CRTmatrix failed" << endl;
    //assert(T2==T3);

    cout << "!"<< endl;
    printAllTimers(cout);

#ifdef CodeBlocks
    cin.get();
#endif
}


