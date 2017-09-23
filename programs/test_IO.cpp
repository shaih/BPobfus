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
/* Test_IO: Test each read and write for the different classes
 */
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <stdexcept>
#include <sys/stat.h>
#include <unistd.h>
#include "../obfus.h"
using namespace std;

#include "../GGH15.h"
#include "../TDMatrixParams.h"


void test_IO(TDMatrixParams& params, long N)
{

    //read and write TDmatrixParams and compare
    std::string dirName = "testIO";

     system(("rm -rf "+dirName).c_str());

#if defined(_WIN32) // Windows mkdir takes only one parameter
    int ret = mkdir(dirName.c_str());
#else               // *nix mkdir takes two parameters
    int ret = mkdir(dirName.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#endif

    ret = chdir(dirName.c_str());

    std::string filename = "testIOFile";

    long bSuccess = true;

    FILE* handle = fopen(filename.c_str(), "wb");
    params.writeToFile(handle);
    fclose(handle);
    handle = fopen(filename.c_str(), "rb");
    TDMatrixParams params2;
    params2.readFromFile(handle);
    fclose(handle);
    chdir ("..");

    if (!(params==params2))
    {
        cout << "error write/read TDMatrixParams class ";
        bSuccess = false;
    }
//    assert(params==params2);

    dirName = "test2";

    //test the function in obfus
    saveParamstoNewDir(dirName, params);
    ret = chdir(dirName.c_str());
    handle = fopen("params.dat", "rb");
    TDMatrixParams params3;
    params3.readFromFile(handle);
    fclose(handle);
    chdir ("..");
    //assert(params==params3);

    if (!(params==params3))
    {
        cout << "error write/read TDMatrixParams class ";
        bSuccess = false;
    }
    dirName = "testIO";
    ret = chdir(dirName.c_str());

//test CRTMatrix writing
    CRTmatrix mCRT(&params);
    mCRT.randomFill(params.mBar,params.wLen);
    handle = fopen(filename.c_str(), "wb");

    mCRT.writeToFile(handle);
    fclose(handle);
    handle = fopen(filename.c_str(), "rb");

    CRTmatrix mCRT2(&params);
    mCRT2.readFromFile(handle, &params);
    fclose(handle);
    ret = chdir ("..");
    //assert(mCRT==mCRT2);

    if (!(mCRT==mCRT2))
    {
        cout << "error write/read CRTmatrix class ";
    }

    //test TDMatrix writing
    TDMatrix tdmp(params);
    ret = chdir(dirName.c_str());
    handle = fopen(filename.c_str(), "wb");
    tdmp.writeToFile(handle);
    fclose(handle);
    handle = fopen(filename.c_str(), "rb");
    TDMatrix tdmp2(params);
    tdmp2.readFromFile(handle, &params);
    fclose(handle);
    ret = chdir ("..");
    if (!(tdmp==tdmp2))
    {
        cout << "error write/read TDMatrix class ";
        bSuccess = false;
    }

    mat_l m1, m2;

    m1.SetDims(params.mBar,params.wLen);
    for (long row = 0; row < params.mBar; row++)
        for (long col = 0; col < params.wLen; col++)
            m1[row][col] = rand() % 100;

    ret = chdir(dirName.c_str());
    handle = fopen(filename.c_str(), "wb");
    writeToFile(m1,handle);
    fclose(handle);
    handle = fopen(filename.c_str(), "rb");
    readFromFile(m2,handle);
    fclose(handle);
    ret = chdir ("..");
    if (!(m1==m2))
    {
        cout << "error write/read mat_l class ";
        bSuccess = false;
    }



    GGH15mmaps maps(params, N); // an (N+1)-node instance
    GGH15mmaps maps2(params);

    ret = chdir(dirName.c_str());
    handle = fopen(filename.c_str(), "wb");
    maps.writeToFile(handle);
    fclose(handle);
    handle = fopen(filename.c_str(), "rb");
    maps2.readFromFile(handle,&params);
    fclose(handle);
    ret = chdir ("..");
    if ((maps!=maps2))
    {
        cout << "error write/read GGH15mmaps class ";
        bSuccess = false;
    }

    if (bSuccess==true)
    cout << "I/O test successful!" << endl;
    else
       cout << "I/O test failed!" << endl;

 }

/************************** main() program below ********************/
/********************************************************************/
#include "../utils/argmap.h"

//#define NFACTORS 8
//static long tinyFactors[NFACTORS] = { 2, 3, 5, 7, 11, 13, 17, 19 };

int main(int argc, char *argv[])
{
#ifdef DEBUG
    NTL::SetSeed((unsigned char*)"test_GGH15",10);
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
    long NN=3;
    amap.arg("NN", NN, "the multi-linearity degree");

    amap.parse(argc, argv); // parses and overrides initail values
    // for each parameter, one per line,

    // assert (k <= NFACTORS); // for this test we don't want many factors
    assert (k <= TDMATRIX_NUM_SMALL_FACTORS);

    if (m < n*(2+k*e)) m = n*(2+k*e);

    TDMatrixParams params(n,k,e,m);

    cout << "NN=" << NN << ", n=" << n << ", m=" << m << ", k="
         << k     << ", e=" << e << ", q=" << params.getQ() << endl;

    test_IO(params,NN);

#ifdef CodeBlocks
    cin.get();
#endif
}
