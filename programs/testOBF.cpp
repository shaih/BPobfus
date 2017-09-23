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
/********************************************************************
test_OBF.cpp - test program for the obfus function implementation

usage Pattern:

long ret=testBP(dirName, trans);

trans is a 2x2 matrix of matrices.
dirname is the directory to which the obfuscation matrices will be written to.

or: run testOBF with command line arguments: dim=#, sig=#, L = #
Where the filename for the transitions is BPs/dim#_sig#_L#.txt
 ************************************************************************/
#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <sys/stat.h>
#include <unistd.h>

#if defined(__unix__) || defined(__unix) || defined(unix)
#include <sys/time.h>
#include <sys/resource.h>
#endif

#include <NTL/BasicThreadPool.h>

#include "../GGH15.h"
#include "../obfus.h"
#include "../utils/argmap.h"

static bool testBP(const std::string& dirName,
		   const Mat< Mat<long> >& trans,
		   long sec, bool pipe, long nThreads=1);


/************************** main() program below ********************/
/********************************************************************/



int main(int argc, char *argv[])
{
    ArgMapping amap;

    std::string name = "P";
    amap.arg("name", name, "name of BP");
    long dim=2;
    amap.arg("dim", dim, "dimension of transition matrices");
    long sig=2;
    amap.arg("sig", sig, "size of the alphabet");
    long L=3; //3
    amap.arg("L", L, "length of BP");
    long sec=80; //3
    amap.arg("sec", sec, "security parameter");
    long nTrials= 1;
    amap.arg("nTrials", nTrials, "number of trials to run");
    std::string dirName = "test";
    amap.arg("testDir", dirName, "where to write the obfuscated program");
    bool pipe = false;
    amap.arg("pipe", pipe, "pipe-based I/O");
    long threads = 1;
    amap.arg("threads", threads, "how many threads to run");


    amap.parse(argc, argv); // parses and overrides initail values

    cout << "\ndim="<<dim<<", sig="<<sig<<", L="<<L<<", sec="<<sec<<endl;

    std::string fullFilename = "BPs/" + name + "_dim" + ToString(dim)
                               + "_sig" + ToString(sig)
                               + "_L" + ToString(L) + ".txt";

    std::fstream fs;
    fs.open(fullFilename.c_str()); // E.g., "BPs/P_dim2_sig2_L3.txt"
    if (!fs.is_open()) // try to create it
    {
      makeRandomBP(dirName, name, dim, sig, L);
      fs.open(fullFilename.c_str());
    }
    if (!fs.is_open()) // failed
    {
        std::cout << "Cannot open input file "<< fullFilename << endl;
#ifdef CodeBlocks
        cin.get();
#endif
        exit(0);
    }

    Mat< Mat<long> > trans; // read the transition matrices from file
    fs >> trans;
    fs.close();

    // Setup a htread pool with almost as many threads as there are cores
    if (threads<1) threads=1;
    cout << "Using "<<threads<< " threads\n";


    bool ret;
    for (long i=0; i<nTrials; i++)
    {
        cout << "started test#" << i << endl;

        if (!(ret=testBP(dirName, trans, sec, pipe, threads)))
        {
            cout << "**test #"<<i<<" failed!!\n";
            break;
        }
    }
    if (ret) cout << "**All "<<nTrials<<" tests succeeded!\n";

#ifdef CodeBlocks
    cin.get();
#endif
}

static bool testBP(const std::string& dirName, const Mat< Mat<long> >& trans,
		   long sec, bool pipe, long nThreads)
{
    system(("rm -rf "+dirName).c_str());

    resetAllTimers();
    initObfuscateNFA(dirName, trans, sec, NULL, pipe, nThreads);

    cout << endl;
    printAllTimers(cout);
#if (defined(__unix__) || defined(__unix) || defined(unix))
    struct rusage rusage;
    getrusage( RUSAGE_SELF, &rusage );
    cout << "  rusage.ru_maxrss="<<rusage.ru_maxrss << endl;
#endif

    resetAllTimers();
    cout << "evaluating" << endl << flush;

    long L = trans.NumRows(); //number of edges

    long two2n = 1L << L;     //=2^L, # of possible strings
    if (two2n>32) two2n=32;   // don't try too many strings

    Vec<long> bitString(INIT_SIZE, L);
    for (long ls=0; ls<two2n; ls++)
    {
        long bit = ls & 1; // lsb
        bitString[0] = bit;
        string ss = ToString(bit);

        mat_l tMat = trans[0][bit];
        for (long i=1; i<L; i++)
        {
            bit = (ls >> i) & 1; // next bit
            bitString[i] = bit;
            ss +=  ToString(bit);

            tMat *= trans[i][bit];
        }

        bool b1 = isZero(tMat);
        bool b4 = evalNFA(dirName, bitString, NULL, nThreads);
        assert(b1==b4);
    }
    printAllTimers(cout);
#if (defined(__unix__) || defined(__unix) || defined(unix))
    getrusage( RUSAGE_SELF, &rusage );
    cout << "  rusage.ru_maxrss="<<rusage.ru_maxrss << endl;
#endif
    return true;
}
