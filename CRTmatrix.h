#ifndef _CRTMATRIX_H_
#define _CRTMATRIX_H_
/* Copyright (C) 2017 IBM Corp.
 *  Licensed under the Apache License, Version 2.0 (the "License"); 
 * you may not use this file except in compliance with the License. 
 * You may obtain a copy of the License at
 *     http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, 
 * software distributed under the License is distributed on an
 * "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
 * either express or implied. See the License for the specific
 * language governing permissions and limitations under the License. 
 */
/*****************************************************************************
 * CRTmatrix.h - A matrix in CRT representation

  Usage Pattern:

  TDMatrixParams params(...); // initialize somehow
  CRTMatrix crtM(params);     //initializes the CRT matrix,
                             // using the factors in params->kFactors

  Converting to and from NTL::mat_ZZ matrices:

  // reconstruct the matrix in ZZ fromat, returning in q the modulus
  NTL::mat_ZZ mTo;
  ZZ q;
  convert(mTo, NTL::ZZ& q, const CRTmatrix& from);

  // reconstruct the matrix in CRT fromat
  void convert(CRTmatrix& to, const mat_ZZ& from, const TDMatrixParams& prms);
 *****************************************************************************/
#include <cstdio>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_lzz_p.h>
#include "mat_l.h"


class TDMatrixParams; // forward decleration

typedef NTL::Vec<NTL::vec_zz_p> CRTvector;
/*------------------------------------------------------------------*/
// A matrix in CRT representaion
class CRTmatrix: public NTL::Vec<NTL::mat_zz_p>
{
public:
    TDMatrixParams *params;

#if 0
    CRTmatrix(const CRTmatrix& other); // copy constructor
#endif
    CRTmatrix(TDMatrixParams* p=NULL)
    {
        params=p;
    }

    // Choose a random n-by-m matrix in CRT representation
    void randomFill(long n, long m);

    long NumRows() const
    {
        return ((length()>0)? (*this)[0].NumRows() : 0);
    }

    long NumCols() const
    {
        return ((length()>0)? (*this)[0].NumCols() : 0);
    }

    // Get a CRT column: crtCol[i] is the j'th column of this[i]
    void getColumn(CRTvector& crtCol, long j);
    CRTvector getColumn(long j)
    {
        CRTvector c;
        getColumn(c, j);
        return c;
    }

    // Convert Mat<long> to CRT format, we assume that params are already set
    CRTmatrix& operator=(const mat_l& M);         // A := M
    CRTmatrix& operator=(const CRTmatrix& other); // A := M

    CRTmatrix& leftMultBy(const vec_l& u);     // A := u*A
    CRTmatrix& leftMultBy(const mat_l& M);     // A := M*A
    CRTmatrix& leftMultBy(const CRTmatrix& M); // A := M*A

    CRTmatrix& operator*=(const CRTmatrix& M); // A := A*M
    CRTmatrix& operator*=(const mat_l& M);     // A := A*M

    CRTmatrix& operator+=(const CRTmatrix& M); // A := A+M
    CRTmatrix& operator+=(const mat_l& M);     // A := A+M

    CRTmatrix& operator-=(const CRTmatrix& M); // A := A-M
    CRTmatrix& operator-=(const mat_l& M);     // A := A-M

    // Invert M modulo the product of the factors, return false on failure
    bool invert(const mat_l& M);

    //invert A, raises error on failure
    bool invert(const CRTmatrix& A);

    // Set as the identity matrix
    void identity(long n);

    // Check that the entries are all "much smaller" than q, where "much smaller"
    // means smaller by a factor of at least 2^bitGap (in absolute value)
    bool isSmall(long bitGap=10) const;

    // binary I/O, returns # of bytes read
    long writeToFile(FILE* handle) const;
    long readFromFile(FILE* handle, TDMatrixParams* prmBuf=NULL);
};

bool operator==(const CRTmatrix& A, const CRTmatrix& B);
inline bool operator!=(const CRTmatrix& A, const CRTmatrix& B)
{
    return !(A==B);
}

// ostream& operator<< (ostream &s, const CRTmatrix& A);

// Choose a random n-by-m matrix in CRT representation
inline void
randomFill(CRTmatrix& M, long n, long m, const TDMatrixParams& prms)
{
    M.params = (TDMatrixParams*) &prms;
    M.randomFill(n,m);
}

// reconstruct the matrix in ZZ fromat, returning in q the modulus
void convert(NTL::mat_ZZ& to, NTL::ZZ& q, const CRTmatrix& from);

// reconstruct the matrix in CRT fromat
void convert(CRTmatrix& to, const mat_ZZ& from, const TDMatrixParams& prms);

inline void conv(CRTmatrix& A, const mat_l& M)
{
    A = M;
}

void generateMatrixPair(CRTmatrix& A, CRTmatrix& Ainv, long n);
void generateMultiPair(CRTmatrix& A, CRTmatrix& Ainv, long n); //function for timing testing purposes

/********************************************************************/
// A conveneince class of a matrix tagged by a number, can be used
// to store the A matrix of GGH15 nodes. Value of tag<0 indicates a
// non-initialized object

class TaggedCRTmatrix
{
public:
    long tag;
    CRTmatrix data;

    TaggedCRTmatrix(TDMatrixParams* p=NULL): data(p)
    {
        tag=-1;
    }

    bool initalized() const
    {
        return (tag >= 0);
    }
};



#endif // _CRTMATRIX_H_
