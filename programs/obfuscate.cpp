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
Usage: ./initialize_x [ name=value ]...
	dim 	dimension of transition matrices  [ default=2 ]
	L 	length of BP  [ default=3 ]
	e 	Power of factors  [ default=3 ]
	sec 	security parameter  [ default=80 ]
	pipe 	pipe-based output  [ default=1 ]
	threads number of threads  [ default=#cores-2 ]
	testDir where to write the obfuscation parameters  [ default=test ]

Uses the 'params.dat' file to determines the parameters, and reads
the secret encoding key from the testDir directory, as well as the
program to be obfuscated, and writes the obfuscated program to the
same directory.
*************************************************************************/
#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <NTL/matrix.h>

#if defined(__unix__) || defined(__unix) || defined(unix)
#include <sys/time.h>
#include <sys/resource.h>
#endif

#include "../utils/argmap.h"
#include "../GGH15.h"
#include "../obfus.h"

//#define DEBUG

NTL_CLIENT
/********************************************************************/

int main(int argc, char *argv[])
{
    ArgMapping amap;

     bool pipe = true;
    long threads = 0;


    std::string name = "P";
    amap.arg("name", name, "name of BP");
    long dim=2;
    amap.arg("dim", dim, "dimension of transition matrices");
    long L=3;
    amap.arg("L", L, "length of BP");
    long sig= 2;
    amap.arg("sig", sig, "alphabet size");
    amap.arg("pipe", pipe, "pipe-based output");
    amap.arg("threads", threads, "number of threads", "#cores-2");

    if (threads <= 0)
        threads = 1;

    long sec = 80;
    std::string dirName = "test";
    amap.arg("testDir", dirName, "where to write the obfuscated program");
#ifdef DEBUG
    sec = 20;
#endif
    amap.arg("sec", sec, "security parameter");

    amap.parse(argc, argv); // parses and overrides initail values

    // Read the transitions matrices from file

    std::string filename = dirName + "/" + name + "_dim" + ToString(dim)
                           + "_sig" + ToString(sig)
                           + "_L" + ToString(L) + ".txt";

    std::fstream fs;
    fs.open(filename.c_str()); // E.g., "BPs/P_dim2_L3_sym2.txt"
    if (!fs.is_open()) // try to create it
    {
        std::cout << "BP file does not exist, creating a random BP in "<< filename << endl;
        //makeRandomBP(dirName, name, dim, sig, L);
        makeRandomBP(dirName, name, dim, sig, L);
        //makeRandomBP(name, dim, sig, L);
        fs.open(filename.c_str());
    }

    if (!fs.is_open())
    {
        std::cout << "Cannot open input file "<< filename << endl;
        return 0;
    }
    Mat< Mat<long> > trans(INIT_SIZE, L, sig);
    fs >> trans;

    cout << filename << " read for obfuscation" << endl;

    //std::string dirName = "mmaps_dim" + ToString(dim) + "_L" + ToString(L) + "_sec" + ToString(sec);

    //read params

    std::string prmsFile = dirName + "/params.dat";
    FILE* handle = fopen(prmsFile.c_str(), "rb");

    if (handle==NULL)
    {
        cout << "missing file: " << prmsFile << ", please run initialize_x first, exiting program" << endl;
        return 0;
    }


    TDMatrixParams p;
    p.readFromFile(handle);
    fclose(handle);

    cout << "m=" << p.m << ", kFactors=" << p.kFactors
	 << ", n=" << p.n << ", e=" << p.e << ", dim=" << dim << endl;

    Vec<double> obfusTimers;

    // Do the actual obfuscation
    obfuscateNFASavedNodes(dirName, trans, &p, pipe, threads);

        const FHEtimer* obfusTimer = getTimerByName("obfuscateNFASavedNodes");
    obfusTimers.append(obfusTimer->getTime());



    printAllTimers(cout);
#if (defined(__unix__) || defined(__unix) || defined(unix))
    struct rusage rusage;
    getrusage( RUSAGE_SELF, &rusage );
    cout << "  rusage.ru_maxrss="<<rusage.ru_maxrss << endl;
#endif

        cout << "--------------" << endl;
    cout << "obfuscating =" << obfusTimers << endl;
    cout << "--------------" << endl;

#ifdef CodeBlocks
    cin.get();
#endif

    return true;
}


