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
/************************************************************************
Usage: ./testOBF_x [ name=value ]...
	name 	name of BP  [ default=P ]
	dim 	dimension of transition matrices  [ default=2 ]
	sig 	size of the alphabet  [ default=2 ]
	L 	length of BP  [ default=3 ]
	sec 	security parameter  [ default=80 ]
	nTrials number of trials to run  [ default=1 ]
	testDir where to write the obfuscated program  [ default=test ]
	pipe 	pipe-based I/O  [ default=0 ]
	threads how many threads to run  [ default=1 ]
*************************************************************************/
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <chrono>
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
#include "../pipe.h"
#include "../utils/argmap.h"

NTL_CLIENT

//#define DEBUG
//#define PRINTDOT

int main(int argc, char *argv[])
{
    long dim=2;
    long L=3;
#ifdef DEBUG
    long sec = 0;
#else
    long sec = 80;
#endif
    bool pipe = true;
    long threads = 0;
    long e=3;

    ArgMapping amap;
    amap.arg("dim", dim, "dimension of transition matrices");
    amap.arg("L", L, "length of BP");
    amap.arg("e", e, "Power of factors");
    amap.arg("sec", sec, "security parameter");
    amap.arg("pipe", pipe, "pipe-based output");
    amap.arg("threads", threads, "number of threads", "#cores-2");
    std::string dirName = "test";
    amap.arg("testDir", dirName, "where to write the obfuscation parameters");
    amap.parse(argc, argv); // parses and overrides initail values
    if (threads<1) threads = numCores()-2;
    if (threads<1) threads = 1;

    //std::string dirName = "mmaps_dim" + ToString(dim) + "_L" + ToString(L)
    //  + "_sec" + ToString(sec);

    system(("rm -rf "+dirName).c_str());
    resetAllTimers();

    TDMatrixParams p;
    GGH15basicParams bp(COMP_GGH15PRMS, dim, L, sec, e);
    p.init(bp.n, bp.k, e, bp.m);//bp.e, bp.m);

    cout << "m=" << p.m << ", kFactors=" << p.kFactors
	 << ", n=" << p.n << ", e=" << p.e << ",dim=" << dim << endl;

  Vec<double> initTimers;

    {FHE_NTIMER_START(totalInitTime);

    cout << "Initializing" << endl;
    bool bsuccess=false;
    if (pipe) {
      bsuccess = pipedInit(dirName, p, 2*L, threads);
    } else {
      SetNumThreads(threads);
      cerr << "using "<<threads<<" threads\n";
      bsuccess = initAllNodes(dirName, p, 2*L);
    }
    if (bsuccess)
      cout << "Initialization successful\n";
    else
      cout << "Initialization failed\n";
    }

    printAllTimers(cout);

        const FHEtimer* initTimer = getTimerByName("totalInitTime");

    initTimers.append(initTimer->getTime());

        cout << "--------------" << endl;
    cout << "Initializing time=" << initTimers << endl;
    cout << "--------------" << endl;


#if (defined(__unix__) || defined(__unix) || defined(unix))
    struct rusage rusage;
    getrusage( RUSAGE_SELF, &rusage );
    cout << "  rusage.ru_maxrss="<<rusage.ru_maxrss << endl;
#endif
    return true;

#ifdef CodeBlocks
    cin.get();
#endif
}
