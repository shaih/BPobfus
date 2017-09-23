#ifndef _GGH15_H_
#define _GGH15_H_
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
/* GGH15.h - an implementation of the graph-based multilinear maps from
 *   [Gentry-Gorbunov-Halevi, TCC 2015].
 *
 * Usage pattern:
 *
 * long N;        // number of nodes in graph
 * TDMatrixParams params(...); // initialize somehow
 *
 * GGH15mmaps maps(params, N); // an N-node instance
 *
 * Mat< mat_l > S;
 * // These are the "plaintext matrices" to be encoded, the row S[i]
 * // will holds the matrices that should be encoded relative to the
 * // edge i -> i+1. Say that we have at least two matrices in each S[i].
 *
 * // Set some values for the S matrices: Each S[i][j] is a matrix with
 * // small entries of dimension n-by-n matrix (where n = params.n).
 * // (Note that n and N are two different numbers.)
 *
 * // Encode each matrix S[i][j] via the encoding C[i][j]
 * Mat< GGH15encoding > C(INIT_SIZE, S.NumRows(), S.NumCols());
 * for (long i=0; i<S.NumRows(); i++) for (long j=0; j<S.NumCols(); j++) {
 *    maps.encode(C[i][j], S[i][j], i, i+1)
 *    // encode the matrices in S, the encoding matrices are in C
 * }
 *
 * // Compute some operations on the encoded matrices
 *
 * for (long i=0; i<N; i++)
 *     maps.subEnc(C[i][0], C[i][1]); // C[i][0] -= C[i][1]
 *     // can add/subtract encoding wrt the same edge
 *
 * for (long i=N-1; i>0; i--) {
 *     maps.mulEnc(C[i-1][0], C[i][0]); // C[i-1][0] *= C[i][0]
 *     // can mult encoding wrt i->j by one wrt j->k,
 *     // the result is an encoding wrt i->k
 *
 *     // After each iteration, C[i-1][0] is encoding wrt (i-1)->N
 * }
 *
 * // Can test for zero encodings wrt 0 -> N
 *
 * if (maps.zeroTest(C[0][0])) std::cout << "Result is zero" << endl;
 * else                        std::cout << "Result is nonzero" << endl;
 */
#include <vector>
#include <cassert>
//#include <mpfr.h>

#include "mat_l.h"         // single-percision integer matrices
#include "pipe.h"

NTL_CLIENT
#include "DGaussSampler.h" // discrete Gaussian sampler
#include "TDMatrix.h"      // matrices with trapdoors

extern bool bGGH15threading;


class  GGH15basicParams  {
public:
  class COMP_PARAMS_TYPE {}; // an indicator type

  long n, e, k, m;

  GGH15basicParams(long nn, long kk, long ee=3, long mm=0) // manually set params
  { n=nn; k=kk; e=ee; m=mm; }

  // Compute params from higher-level constraints
  GGH15basicParams(COMP_PARAMS_TYPE comp, long nStates, long diam,
       	           long secParam=80, long ee=3);
};
// An indicator variable, to be used in the constructor
const GGH15basicParams::COMP_PARAMS_TYPE COMP_GGH15PRMS
      = GGH15basicParams::COMP_PARAMS_TYPE();

/********************************************************************/
// One node in the GGH15 encoding graph, contain an n-by-m matrix A with its
// trapdoor, as well as two random transformation matrices T and T^{-1} (mod q)
class GGH15node {
  long index;
  TDMatrixParams* params; // Parameters (e.g. CRT basis)
  TDMatrix A;
  mat_l P;           // inner transformation matrices, P is small
  CRTmatrix Pinv;    //                                but P^{-1} is not
  CRTmatrix T, Tinv; // outter transformation matrices in CRT representation

public:
  // Empty constructor, makes it easier to allocate a vector of these
   //GGH15node(): params(NULL) {}
   GGH15node(TDMatrixParams* p=NULL) : params(p), A(p), Pinv(p), T(p), Tinv(p) {}

  // Constructor: Choose random A w/ trapdoor, and random T, Tinv = T^{-1}
   GGH15node(TDMatrixParams& p, long idx) { initGGH15node(p, idx); }

  void initGGH15node(TDMatrixParams& p, long idx);

  // Access methods
  const TDMatrixParams* getPrms() const { return params; }

  const TDMatrix& getA() const { return A; }
  void getCRTA(CRTmatrix& crtA) const { A.getA(crtA); }
  long getIndex() const {return index;}

  const mat_l& getP() const { return P; }
  const CRTmatrix& getPinv() const { return Pinv; }

  const CRTmatrix& getT() const { return T; }
  const CRTmatrix& getTinv() const { return Tinv; }

  // A wrapper, still unused
  long writeToFile(FILE* handle) {return writeToActualFile(handle);}

  long writeToActualFile(FILE* handle); // The actual thing

  long readFromFile(FILE* handle, TDMatrixParams* prmBuf);
};
/********************************************************************/

  bool operator==(const GGH15node& A, const GGH15node& B);
  inline bool operator!=(const GGH15node& A, const GGH15node& B){return (!(A==B));}

/********************************************************************/
// A GGH15-encoded matrix
class GGH15encoding {
  friend class GGH15mmaps;

  long from, to;
  CRTmatrix data;

public:
  void leftMultBy(const GGH15encoding& other) {
    assert(from==other.to);
    data.leftMultBy(other.data);
    from=other.from;
  }

  GGH15encoding& operator*=(const GGH15encoding& other) {
    assert(to==other.from);
    data *= other.data;
    to=other.to;
    return *this;
  }
  GGH15encoding& operator+=(const GGH15encoding& other) {
    assert(to==other.to && from==other.from);
    data += other.data;
    return *this;
  }
  GGH15encoding& operator-=(const GGH15encoding& other) {
    assert(to==other.to && from==other.from);
    data -= other.data;
    return *this;
  }

  // access methods, used only for debugging
  long fromNode() const { return from; }
  long toNode() const { return to; }
    void setFrom(long Fr)  { from=Fr;}
  void setTo(long T)  { to=T; }
   CRTmatrix& getData()  { return data; }

  // binary I/O, returns # of bytes read/written
  long writeToFile(FILE* handle);
  long readFromFile(FILE* handle, TDMatrixParams& prmBuf);
};
/********************************************************************/


/********************************************************************/
// An instance of multilinear maps, with all its secret parameters
class GGH15mmaps {
  TDMatrixParams *params;

  // The GGH15 secret params include the trapdoor, T,T^{-1}, and P,P^{-1}
  NTL::Vec< GGH15node > nodes; // Nodes 1..n-1 with their trapdoors

  CRTmatrix An;// The An matrix (actually n-by-1 "matrix")

 public:
  // initialize a GGH15 instance with n vertices
  GGH15mmaps(TDMatrixParams& p, long n);
  //GGH15mmaps(TDMatrixParams* p=NULL) : params(p), An(p) {}
  GGH15mmaps(TDMatrixParams& p) : params(&p), An(&p) {}
  // We only keep nodes 1..n-1, nodes 0 and n are not real nodes.
  long numNodes() const { return nodes.length()+2; }

  // Encode plaintext matrix wrt the path i -> j. The toA argument is
  // an optimization:
  //  + If toA is non-NULL and toA->tag==j, then it is assumed that
  //    toA->data is the matrx A_j.
  //  + If toA is non-NULL and toA->tag!=j, then toA->data is set to
  //    A_j and toA->tag is set to j.
  bool encode(GGH15encoding& C, const mat_l& S, long i, long j,
              TaggedCRTmatrix* toA=NULL);

  // two- and three-argument variants
  void addEnc(GGH15encoding& arg1, const GGH15encoding& arg2)
  { arg1 += arg2; }
  void addEnc(GGH15encoding& target,
	      const GGH15encoding& arg1, const GGH15encoding& arg2)
  { target = arg1; target += arg2; }

  void subEnc(GGH15encoding& arg1, const GGH15encoding& arg2)
  { arg1 -= arg2; }
  void subEnc(GGH15encoding& target,
	      const GGH15encoding& arg1, const GGH15encoding& arg2)
  { target = arg1; target -= arg2; }

  // two- and three-argument variants
  void mulEnc(GGH15encoding& arg1, const GGH15encoding& arg2)
  { arg1 *= arg2; }
  void mulEnc(GGH15encoding& target,
	      const GGH15encoding& arg1, const GGH15encoding& arg2)
  { target = arg1; target *= arg2; }

  bool zeroTest(const GGH15encoding& enc)
  { return (enc.from==0 && enc.to==numNodes()-1 && enc.data.isSmall()); }

  TDMatrixParams* getParams() const { return params; }
  const GGH15node& getNode(long index) const {return (nodes[index]);}
  const CRTmatrix& getAn() const {return An;}

  long writeToFile(FILE* handle);
  long readFromFile(FILE* handle, TDMatrixParams* prmBuf);

  // A debugging method to check that an encoding was done right
  bool verifyEncoding(const GGH15encoding& C, const mat_l& S) const;

};
/********************************************************************/

  bool operator==(const GGH15mmaps& A, const GGH15mmaps& B);
  inline bool operator!=(const GGH15mmaps& A, const GGH15mmaps& B){return (!(A==B));}

/********************************************************************/
class GGH15ptxt  {
  long from,to;
  std::string tag;
  mat_l data;

public:
 GGH15ptxt(): tag(), data() { from=to=0; }

  GGH15ptxt(long f, long t, const std::string& tstr, const mat_l& d)
  { from = f; to = t; data = d; tag = tstr; }

  long getTo() const    { return to; }
  long getFrom() const  { return from; }
  const std::string& getTag() const { return tag; }
  const mat_l& getData() const { return data; }
  mat_l& getData() { return data; } // allows non-const matrix

  void setTo(long t)    { to=t; }
  void setFrom(long f)  { from=f; }
  void setTag(const std::string& s) { tag=s; }
  void setData(const mat_l& m) { data=m; }
};


typedef NTL::Pair<FILE*,GGH15node*> nodePair;
bool initAllNodes(std::string dirName, TDMatrixParams& p,
		  long numNodes, SynchronizedPipe<nodePair>* pipe=NULL);
//initialize all nodes without saving to an mmaps structure


bool encodeMatrix(const std::string& dirName,
		  GGH15encoding& C, const mat_l& S, long i, long j,
                  GGH15node* iNode, GGH15node* jNode, CRTmatrix& An,
                  TaggedCRTmatrix* toAptr, TDMatrixParams *params, int nNodes);

#endif // ifndef _GGH15_H_
