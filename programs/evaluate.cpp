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
/* evaluate.cpp - evaluate an obfuscated BP on a single input
 *
 * Usage: ./evaluate_x [ name=value ]...
 *	name 	name of BP  [ default=P ]
 *	dim 	dimension of transition matrices  [ default=2 ]
 *	L 	langth of BP  [ default=3 ]
 *	threads number of threads  [ default=#cores-2 ]
 *	sig 	alphabet size  [ default=2 ]
 *	testDir where to write the obfuscated program  [ default=test ]
 *	input 	input string of size L  [ default=all zero ]
 *
 * The program reads the obfuscated program from the directory testDir
 * and the input from the command-line parameters, and performs the
 * evaluation.
 *************************************************************/
#include <fstream>
#include <iostream>
#include <string>

#if defined(__unix__) || defined(__unix) || defined(unix)
#include <sys/time.h>
#include <sys/resource.h>
#endif

#include <NTL/matrix.h>

#include "../utils/argmap.h"
#include "../GGH15.h"
#include "../obfus.h"

NTL_CLIENT
/********************************************************************/

static bool verifyInputString(const std::string& input, long L, long nSym);

int main(int argc, char *argv[])
{
    ArgMapping amap;
    long threads = 0;


    std::string name = "P";
    amap.arg("name", name, "name of BP");
    long dim=2;
    amap.arg("dim", dim, "dimension of transition matrices");
    long L=3;
    amap.arg("L", L, "langth of BP");
    amap.arg("threads", threads, "number of threads", "#cores-2");
    long sig= 2;
    amap.arg("sig", sig, "alphabet size");
    std::string dirName = "test";
    amap.arg("testDir", dirName, "where to write the obfuscated program");

    std::string input="";
    amap.arg("input", input, "input string of size L", "all zero");

    amap.parse(argc, argv); // parses and overrides initail values

    if (input.length()!=L)
    {
        input = "";
        for (long i = 0; i<L; i++)
            input+="0";
    }

    cout << "L="<<L<<", sig="<<sig<<", dim=" << dim << endl;

    // verify that the input string is valid
    if (!verifyInputString(input, L, sig))
    {
        std::cerr << "invalid input string "<<input
                  << " for L="<<L<<", sig="<<sig<<", dim=" << dim << endl;

        return 0;
    }

    // from ascii to binary representation
    Vec<long> binInput(INIT_SIZE, L);
    for (long i=0; i<L; i++) binInput[i] = input[i] - '0';

    Vec<double> evalTimers;

    // Evaluate the obfuscated BP
    //std::string dirName = "OBF" + name + "_dim" + ToString(dim)
    //                      + "_sig" + ToString(sig)
     //                     + "_L" + ToString(L);
    //cout << "evaluating" << std::flush;

    bool out = evalNFA(dirName, binInput, NULL, threads);
    cout << name << "("<<input<<")="<<out<<endl;

    /* evaluate 32 options?
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
        for (long i=1; i<L; i++)*
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
    */ //evaluate 32 options

    const FHEtimer* evalTimer = getTimerByName("evalNFA");
    evalTimers.append(evalTimer->getTime());


    // For testing purposes, look for a cleartext version of the BP

    Mat< Mat<long> > trans;
    std::string clearname = dirName + "/" + name + "_dim" + ToString(dim)
                            + "_sig"+ ToString(sig)
                            + "_L" + ToString(L)
                            + ".txt";

   // cout << clearname.c_str() << " = BP file opened" << endl;
    std::fstream clearfile;
    clearfile.open(clearname.c_str()); // E.g., "BPs/P_dim2_L3_sym2.txt"

    if (!clearfile.is_open()) //cleartext is not there, return
    {
#ifdef CodeBlocks
        cin.get();
#endif
        cout << clearname.c_str() << " not found" << endl;
        exit(0);
    }

    if (clearfile.is_open()) // file exists
        clearfile >> trans;           // read cleartext BP

    if (trans.NumRows()>0)    // compare to cleartest evaluation
    {
        long sym = binInput[0]; // first symbol
        mat_l tMat = trans[0][sym];
        for (long i=1; i<L; i++)
        {
            sym = binInput[i];    // next symbol
            tMat *= trans[i][sym];
        }
        bool outClear = isZero(tMat);

        if (out != outClear)
            cerr << "  Error: output should have been "<<outClear<< endl;
    }

    cout << "--------------" << endl;
    cout << "evaluating =" << evalTimers << endl;

    cout << "--------------" << endl;

    cout << "!"<< endl;
    printAllTimers(cout);
#if (defined(__unix__) || defined(__unix) || defined(unix))
    struct rusage rusage;
    getrusage( RUSAGE_SELF, &rusage );
    cout << "  rusage.ru_maxrss="<<rusage.ru_maxrss << endl;
#endif


#ifdef CodeBlocks
    cin.get();
#endif

    return true;
}

static bool verifyInputString(const std::string& input, long L, long sig)
{
    if (L != (long)input.length()) return false; // length should match

    char lastSym = '0' + sig-1; // every symbol must be between 0 and nSym-1
    for (long i=0; i<L; i++)
        if (input[i] < '0' || input[i] > lastSym)
            return false;

    return true;
}
