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
test_timing.cpp - test program for timing the initialization, obfuscation and evaluation

usage Pattern:

run test_timing with command line arguments: dim=#, sig=#, L = #
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

static bool testBP(const std::string& dirName, const Mat< Mat<long> >& trans, long sec,
                   Vec<double> &initTimers, Vec<double> &obfusTimers, Vec<double> &evalTimers, long threads, bool pipe);

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
    bool pipe=0;
    amap.arg("pipe", pipe, "pipe-based I/O");
    std::string dirName = "test";
    amap.arg("testDir", dirName, "where to write the obfuscated program");
    long threads = 1;
    amap.arg("threads", threads, "how many threads to run");

    amap.parse(argc, argv); // parses and overrides initail values

    Vec<double> initTimers;
    Vec<double> obfusTimers;
    Vec<double> evalTimers;

    if (threads<1) threads=1;
    cout << "Using "<<threads<< " threads\n";
    cout << "--------------" << endl;

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

    for (long i=0; i<nTrials; i++)
    {
      if (!testBP(dirName, trans,sec,initTimers, obfusTimers, evalTimers,
		  threads, pipe))
        {
            cout << "**test #"<<i<<" failed!!\n";
            break;
        }
    }

    cout << "--------------" << endl;

    cout << "Initializing =" << initTimers << endl;
    cout << "obfuscating =" << obfusTimers << endl;
    cout << "evaluating =" << evalTimers << endl;

    cout << "--------------" << endl;

#if (defined(__unix__) || defined(__unix) || defined(unix))
    struct rusage rusage;
    getrusage( RUSAGE_SELF, &rusage );
    cout << "  rusage.ru_maxrss="<<rusage.ru_maxrss << endl;
#endif

    std::string prmFile = dirName+"/params.dat";
    FILE* handle = fopen(prmFile.c_str(), "rb");
    if (handle != 0)
    {
        TDMatrixParams p;
        p.readFromFile(handle);
        fclose(handle);
        cout << "m=" << p.m
             << ", kFactors=" << p.kFactors
             << ", n=" << p.n
             << ", e=" << p.e << endl;
    }


#ifdef CodeBlocks
    cin.get();
#endif
}

static bool testBP(const std::string& dirName, const Mat< Mat<long> >& trans,
		   long sec, Vec<double> &initTimers, Vec<double> &obfusTimers,
		   Vec<double> &evalTimers, long threads, bool pipe)
{
    system(("rm -rf "+dirName).c_str());

    resetAllTimers();

#if (defined(__unix__) || defined(__unix) || defined(unix))
      struct rusage rusage;
    getrusage( RUSAGE_SELF, &rusage );
    cout << "  rusage.ru_maxrss="<<rusage.ru_maxrss << endl;
#endif

    initObfuscateNFA(dirName, trans, sec, NULL, pipe, threads);

#if (defined(__unix__) || defined(__unix) || defined(unix))
    getrusage( RUSAGE_SELF, &rusage );
    cout << "  rusage.ru_maxrss="<<rusage.ru_maxrss << endl;
#endif

     cout << "finished obfuscation" << endl;

    const FHEtimer* initTimer = getTimerByName("initAllNodes");//GGH15mmaps");
    initTimers.append(initTimer->getTime());
    const FHEtimer* obfusTimer = getTimerByName("obfuscateNFASavedNodes");
    obfusTimers.append(obfusTimer->getTime());

    printAllTimers(cout);

#if (defined(__unix__) || defined(__unix) || defined(unix))
    getrusage( RUSAGE_SELF, &rusage );
    cout << "  rusage.ru_maxrss="<<rusage.ru_maxrss << endl;
#endif

    resetAllTimers();

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
        bool b4 = evalNFA(dirName, bitString, NULL, threads);
        assert(b1==b4);
    }
    const FHEtimer* evalTimer = getTimerByName("evalNFA");
    evalTimers.append(evalTimer->getTime());

    printAllTimers(cout);

    return true;
}
