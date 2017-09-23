#ifndef _TDMATRIX_H_
#define _TDMATRIX_H_
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
/*
 * TDMatrix.h - Matrixes with trapdoors,
 *   implementing the Efficient algorithm SampleD - Algorithm 3 from MP12.
 *
 * Usage pattern:
 *
 * int res = sampleWithTrapdoor(vec_l &xOut, vec_l& p, vec_l& z,
                                 Vec<vec_zz_p>& U) const


 *   TDMatrixParams params; // initialize parameters somehow
 *   TDMatrix A(params);  // choose a random matrix A with a trapdoor
 *
 *   Vec<double> syndrom;   // some desired syndrom
 *   Vec<long> x;
 *   int res = A.sampleWithTrapdoor(x, syndrom); // returns 0 on suceess
 *       // sample x such that A * x = syndrom (mod q)

 *   Or:

*  int res = A.sampleWithTrapdoor(x, p, z, syndrom); // returns 0 on suceess
 *       // sample x such that A * x = syndrom (mod q)
 *      //this debug version also returns p and z that are chosen randomly
        //s.t. x = p + [R I] ^T * z, where R is the trapdoor chosen randomly as a small matrix

* mat_l TrapdoorR =  getR(); //returns the trapdoor matrix

Initializing the TDMatrix class:

// Choose a (pseudo)random A with a trapdoor R s.t. A*[R/I]=G
// A=(aBar| aPrime ) with aPrime=G-aBar*R isn't computed explicitly here
// R is the trapdoor chosen randomly as a small matrix
TDMatrix tdM;
tdM.initTDmatrix(TDMatrixParams& prms, const CRTmatrix* abar)

The (psuedo) random A matrix used for sampling can be read (for debug purposes) by caling:

CRTmatrix A;
tdM.getA(A);


 */
#include <utility> // define std::pair
#include <NTL/mat_lzz_p.h>

#include "utils/tools.h"
#include "mat_l.h"
#include "TDMatrixParams.h"
#include "DGaussSampler.h"
#include "CRTmatrix.h"

typedef Vec<long> vec_l;

class CRTmatrix; // forward decleration

/*------------------------------------------------------------------*/
/** TDMatrix, the main class for sampling with a trapdoor.
 *
 * The constructor generate a (pseudo)random A and a small R such that
 * A x [R/I] = G (mod q) where [R/I] denote R on top of the identity I.

* Usage Pattern:

* CRTmatrix B; //some matrix
* Vec<vec_zz_p> colMod;
* B.getColumn(colMod, j); // The j'th column of B in CRT format
* vec_l sampledColumn;    // this is where we put the vector x_j
* CRMatrix abar; //another CRMatrix
* TDMatrix A;
* TDMatrixParams prms; //params initiated

* A.initTDmatrix(prms, &abar);
* or
* A.initTDmatrix(prms);

* A.sampleWithTrapdoor(sampledColumn, colMod);// sample a column
* or
* vec_l p, z;
* return sampleWithTrapdoor(xOut, p, z, U);

* this returns the intermediate p and z vectors as well


 **/
class TDMatrix
{
private:
    TDMatrixParams *params; // point to the parameters

    CRTmatrix aBarCRT;  //aBar in chinese remaindering representation
    mat_l trapDoorMatR; //Trapdoor matrix of dimension mBar x wLen

    // A has the form (aBar| aPrime ), where aPrime = G - aBar x R


    // The gaussSamp sampler is used to sample the pertubation vector.
    // the pertubation vector p in [MP'12] is chosen with convariance
    //     SigmaP = sigmaX I - sigmaG [R/I] (R^t|I)
    // where sigmaX, sigmaG are some costants. We then compute from
    // the covariance matrix sigmaP the conditional covarinace, and
    // store them in sigmaPconditional;

    DiscreteGaussianSampler gaussSamp;

public:

    // A static variable for debugging purposes
    // VJS: watch out with threads
    static std::atomic<int> maxSample; // keep the largest sample ever drawn

    //Empty constructor,  makes it easier to initialize a vector of these
    TDMatrix(TDMatrixParams* p=NULL) : params(p), aBarCRT(p) {}

    // Choosing a (pseudo)random A with a trapdoor
    void initTDmatrix(TDMatrixParams& prms, const CRTmatrix* abar=NULL);
    explicit TDMatrix(TDMatrixParams& prms, const CRTmatrix* abar=NULL)
    {
        initTDmatrix(prms, abar);
    }


    void getA(CRTmatrix& A) const; // return the matrix A in CRT format
    const CRTmatrix& getABar() const
    {
        return aBarCRT;
    }
    const mat_l& getR() const
    {
        return trapDoorMatR;
    }
    const TDMatrixParams* getParams() const
    {
        return params;
    }

    // After choosing A,R, we cant to sample x s.t. A*x=u (mod q)

    // A debugging version, returns also p and z
    int sampleWithTrapdoor(vec_l &xOut,
                           vec_l& p, vec_l& z, Vec<vec_zz_p>& U) const;

    int sampleWithTrapdoor(vec_l& xOut, Vec<vec_zz_p>& U) const
    {
        vec_l p, z;
        return sampleWithTrapdoor(xOut, p, z, U);
    }

    long writeToFile(FILE* handle);
    long readFromFile(FILE* handle, TDMatrixParams* prmBuf=NULL);
};

bool operator==(const TDMatrix& A, const TDMatrix& B);
inline bool operator!=(const TDMatrix& A, const TDMatrix& B)
{
    return (!(A==B));
}

// Multiply a vector by the G matrix modulo the current NTL modulus, y=G*x.
// G is and n-by-m matrix with m=n*e*numOfFactors, and is definde by means
// of the vector of factors and the exponent e, each row of G is mostly zero,
// except for a progression
//    ( 1, f1,...,f1^e,  (f1^e f2), (f1^e f2^2), ..., (f1^e...fk^{e-1}) )
template<class T>
void multByG(vec_zz_p &out, const Vec<T> &in,
             long n, const vec_l& factors, long e);
// implemented for Vec<long> and Vec<zz_p>

// A function that returns z such that G*z = v (mod q)
int sampleG(vec_l &xOut, const Vec<vec_zz_p>& v, TDMatrixParams* params);

#endif // _TDMATRIX_H_
