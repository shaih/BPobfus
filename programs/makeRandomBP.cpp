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
/* makeRandomBP.cpp - generate random read-once BPs
 */
#include <fstream>
#include <iostream>
#include <string>

#include <NTL/ZZ.h>
#include <NTL/matrix.h>

#include "../utils/tools.h"
#include "../utils/argmap.h"

NTL_CLIENT

/** Usage example:
 *
 * % ./makeRandomBP_x name=myProgram dim=5 L=4 sig=4
 *
 * This will generate a length-5 BP over 4-symbol alphabet
 * with 5x5 transition matrices
 **/
int main(int argc, char *argv[])
{
    extern bool bGGH15threading;

    ArgMapping amap;

    std::string dirName = "test";
     amap.arg("testDir", dirName, "where to write the obfuscated program");

    std::string name = "P";
    amap.arg("name", name, "Name of BP");
    long dim=2;
    amap.arg("dim", dim, "dimension of transition matrices");
    long L=3;
    amap.arg("L", L, "langth of BP");
    long sig= 2;
    amap.arg("sig", sig, "alphabet size");


    amap.parse(argc, argv); // parses and overrides initail values

    bGGH15threading = 0; //off

    // For each step i=0,1,...,L-1 and each symbol sigma=0,1,...,nSym-1,
    // choose a zero-matrix with probability 1/L, random 0/1 matrix otherwise

    Mat< Mat<long> > trans(INIT_SIZE, L, sig);
    for (long i=0; i<L; i++) for (long iSig=0; iSig<sig; iSig++) {
        Mat<long>& M = trans[i][iSig];
        M.SetDims(dim,dim);
        if (NTL::RandomBnd(L)==0) clear(M);
        else
	  for (long row=0; row<dim; row++) for (long col=0; col<dim; col++) {
              M[row][col] = NTL::RandomBnd(2);
	  }
      }

    // Write the transitions matrices to file

    std::string filename = dirName + "/" + name + "_dim" + ToString(dim)
                                + "_sig" + ToString(sig)
                                + "_L" + ToString(L)
                                + ".txt";

    std::fstream fs;
    fs.open(filename.c_str(), std::ios::out); // E.g., "BPs/P_dim2_sig2_L3.txt"
    if (!fs.is_open())
    {
      std::cout << "Cannot open input file "<< filename << endl;
      exit(0);
    }
    fs << trans << endl;
    fs.close();

#ifdef CodeBlocks
    cin.get();
#endif
}
