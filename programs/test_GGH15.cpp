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
/* test_GGH15.cpp - test program for the GGH15 MMAPs implementation
 */
#include "../GGH15.h"
#include "../obfus.h"
#include "../utils/argmap.h"

#if defined(__unix__) || defined(__unix) || defined(unix)
#include <sys/time.h>
#include <sys/resource.h>
#endif

#include <NTL/BasicThreadPool.h>

static void test1(const Mat< mat_l >& S, TDMatrixParams& params);
static void test2(const Mat< mat_l >& S, TDMatrixParams& params);
static void test3(const Mat< mat_l >& S, TDMatrixParams& params);

void test_GGH15(TDMatrixParams& params, long L)
{
    Mat< mat_l > S(INIT_SIZE, L, 2);
    // These are the "plaintext matrices" to be encoded, the row S[i]
    // will holds the matrices that should be encoded relative to the
    // edge i -> i+1. Say that we have at least two matrices in each S[i].

    // Set some values for the S matrices: Each S[i][j] is a matrix with
    // small entries of dimension n-by-n matrix (where n = params.n).
    for (long i=0; i<L; i++)
    {
        setSmall(S[i][0], params.n, params.n, params.r);
        S[i][1] = S[i][0];
#ifdef DEBUG
        cout << "S["<<i<<"]=" << S[i][0] << endl;
#endif
    }

    test1(S, params);
    test2(S, params);
    clear(S[0][0]);

    clear(S[0][0]);
    test3(S, params);

    cout << "\nmaxSigma="<<Gaussian1Dsampler::maxSigma
         << ", maxSample="<<TDMatrix::maxSample << endl;

    cout << "!"<< endl;
    printAllTimers(cout);

#if (defined(__unix__) || defined(__unix) || defined(unix))
    struct rusage rusage;
    getrusage( RUSAGE_SELF, &rusage );
    cout << "  rusage.ru_maxrss="<<rusage.ru_maxrss << endl;
#endif
}

static void test1(const Mat< mat_l >& S, TDMatrixParams& params)
{
    mat_l zeroMat(INIT_SIZE, params.n, params.n);
    clear(zeroMat);

    long N = S.NumRows();
    Mat< GGH15encoding > C(INIT_SIZE, N, S.NumCols());
    GGH15mmaps maps(params, N+1); // an (N+1)-node instance
    TaggedCRTmatrix toA;
    for (long i=0; i<N; i++) for (long j=0; j<S.NumCols(); j++)
        {
            maps.encode(C[i][j], S[i][j], i, i+1, &toA);
            // encode the matrices in S, the encoding matrices are in C
        }
    GGH15encoding zeroEnc;
    maps.encode(zeroEnc, zeroMat, 1, 2);

    // Also generate an encoding of the matrix Prod_i S[i][0]
    // relative to the edge 0 -> N
    mat_l SS = S[0][0];
    for (long i=1; i<N; i++) SS *= S[i][0];
    GGH15encoding cc;
#ifdef DEBUG
    cout << "\\prod S_i=" << SS << endl;
#endif
    maps.encode(cc, SS, 0, N);

    // Compute some operations on the encoded matrices
    GGH15encoding prod1 = C[0][0];
    GGH15encoding prod2 = C[0][0];
    for (long i=1; i<N; i++)
    {
        // forward multiply: C[0][0] *= C[i][0]
        maps.mulEnc(prod1, C[i][0]);

        if (i==zeroEnc.fromNode())
            maps.mulEnc(prod2, zeroEnc);
        else
            maps.mulEnc(prod2, C[i][0]);

        // backward multiply: C[N-i-1][1] *= C[N-i][1]
        maps.mulEnc(C[N-i-1][1], C[N-i][1]);

        // can mult encoding wrt i->j by one wrt j->k,
        // the result is an encoding wrt i->k
    }

    // Test that the two results encode the same thing, and that they
    // both encode the product of matrices

    maps.subEnc(prod1, cc); // prod1 -= cc
    if (maps.zeroTest(prod1)) std::cout << "PASSED encoding test 1" << endl;
    else                      std::cout << "FAILED encoding test 1" << endl;

    maps.subEnc(cc, C[0][1]);
    if (maps.zeroTest(cc)) std::cout << "PASSED encoding test 2" << endl;
    else                   std::cout << "FAILED encoding test 22" << endl;

    // Check that the product itself is nonzero
    if (maps.zeroTest(prod2)) std::cout << "PASSED encoding test 3" << endl;
    else                      std::cout << "FAILED encoding test 3" << endl;


}

static void test2(const Mat< mat_l >& S, TDMatrixParams& params)
{
    long N = S.NumRows();
    Mat< GGH15encoding > C(INIT_SIZE, N, S.NumCols());
    GGH15mmaps maps(params, 2*N); // an 2N-node instance, two parallel chains

    // encode the same matrices on both chains
    maps.encode(C[0][0], S[0][0], 0, 1);
    maps.encode(C[0][1], S[0][0], 0, 2);
    for (long i=1; i<N-1; i++)
    {
        maps.encode(C[i][0], S[i][0], 2*i-1, 2*i+1);
        maps.encode(C[i][1], S[i][0], 2*i, 2*i +2);
    }
    maps.encode(C[N-1][0], S[N-1][0], 2*N-3, 2*N-1);
    maps.encode(C[N-1][1], S[N-1][0], 2*N-2, 2*N-1);

    GGH15encoding prod1 = C[0][0];
    GGH15encoding prod2 = C[0][1];
    for (long i=1; i<N; i++)
    {
        prod1 *= C[i][0];
        prod2 *= C[i][1];
    }
    prod1 -= prod2;
    if (maps.zeroTest(prod1)) std::cout << "PASSED encoding test 4" << endl;
    else                      std::cout << "FAILED encoding test 4" << endl;

    if (maps.zeroTest(prod2)) std::cout << "FAILED encoding tests 4.5" << endl;
}

static void test3(const Mat< mat_l >& S, TDMatrixParams& params)
{
    params.n += 2;
    params.mBar += 4;
    params.wLen += 2 * params.kFactors * params.e;
    params.m += 4 + 2 * params.kFactors * params.e;

    long N = S.NumRows();
    Mat< GGH15encoding > C(INIT_SIZE, N, S.NumCols());
    GGH15mmaps maps(params, 2*N); // an 2N-node instance, two parallel chains

    GGH15ptxt p1, p2;  // encode randomized matrices on both chains
    for (long i=0; i<N; i++)
    {
        randomizeTransitions(p1, p2, S[i][0], i, 0, N, params.n);
#ifdef DEBUG
        cout << "Encoding for edges "
             << p1.getFrom() <<"->"<< p1.getTo() << ", "
             << p2.getFrom() <<"->"<< p2.getTo() << endl;
        cout << "matrices are:\n" << p1.getData() << endl << p2.getData() << endl;
#endif // DEBUG
        maps.encode(C[i][0], p1.getData(), p1.getFrom(), p1.getTo());
        maps.encode(C[i][1], p2.getData(), p2.getFrom(), p2.getTo());
    }

    GGH15encoding prod1 = C[0][0];
    GGH15encoding prod2 = C[0][1];
    for (long i=1; i<N; i++)
    {
        prod1 *= C[i][0];
        prod2 *= C[i][1];
    }
    prod1 -= prod2;
    if (maps.zeroTest(prod1)) std::cout << "PASSED encoding test 5" << endl;
    else                      std::cout << "FAILED encoding test 5" << endl;
}


/************************** main() program below ********************/
/********************************************************************/


//#define NFACTORS 8
//static long tinyFactors[NFACTORS] = { 2, 3, 5, 7, 11, 13, 17, 19 };

int main(int argc, char *argv[])
{
#ifdef DEBUG
    NTL::SetSeed((unsigned char*)"test_GGH15",10);
#endif
    ArgMapping amap;
    long L=3;
    amap.arg("L", L, "the chain length (in edges)");
    long k=0; //0
    amap.arg("k", k, "number of factors", "0: heuristic");
    long n=5;
    amap.arg("n", n, "the small dimension"); // no default info
    long e=3;
    amap.arg("e", e, "the power of each factor");
    long m=0;
    amap.arg("m", m, "the dimension m", "0: heuristic");
    long threads = 4;
    amap.arg("threads", threads, "how many threads to run");

    amap.parse(argc, argv); // parses and overrides initail values
    // for each parameter, one per line,

    GGH15basicParams bp(COMP_GGH15PRMS, n, L, /*sec=*/80, e);
    if (k>0) bp.k=k;
    if (m > 0) bp.m=m;
    if (m < n*(2+k*e)) m = n*(2+k*e);

    TDMatrixParams params(bp.n, bp.k, e, bp.m);

    cout << "L=" << L << ", n=" << bp.n << ", m=" << bp.m << ", k="
         << bp.k     << ", e=" << e << ", sigmaX=" << params.sigmaX
	 << ", q=" << params.getQ() << endl;

    if (threads<1) threads=1;
    else if (threads > numCores()-2)   // numCores() defined in utils/tools.h
    {
        threads = numCores()-2;
        if (threads<1) threads=1;
    }
    NTL::SetNumThreads(threads);
    cout << "Using "<<threads<< " threads\n";

    test_GGH15(params,L);

#ifdef CodeBlocks
    cin.get();
#endif
}
