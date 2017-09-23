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
/**********************************************************************
 * tools.cpp: misc. utilities used within the main program
 ***********************************************************************/
#include <NTL/version.h>
#include "tools.h"
#include "timing.h"
NTL_CLIENT

// A is an n x n matrix, we attempt to compute its inverse mod p^r.
// If A is invertible, X is set to A^{-1} and return value is true.
// Otherwise, return value is false.
// zz_p::modulus() is assumed to be p^r, for p prime, r >= 1.

bool ppInvert(Mat<zz_p>& X, const Mat<zz_p>& A, long p, long r)
{
#if ((NTL_MAJOR_VERSION>9)||((NTL_MAJOR_VERSION==9) && (NTL_MINOR_VERSION>=7)))
  // for recent enough versions of NTL, just use native calls
  FHE_TIMER_START;
  zz_p d;
  relaxed_inv(d, X, A);    // X = A^{-1}
  return d != 0;

#else
  // for older versions of NTL, use Hensel lifting

  if (r == 1) { // use native inversion from NTL
    zz_p d;
    inv(d, X, A);    // X = A^{-1}
    return d != 0;
  }

  // begin by inverting A modulo p

  // convert to long for a safe transaltion to mod-p objects
  Mat<long> tmp;
  conv(tmp, A);
  {
    zz_pPush push(p);  // temporarily work mod p
    zz_p d;

    Mat<zz_p> A1, Inv1;
    conv(A1, tmp);      // Recover A as a Mat<zz_p> object modulo p
    inv(d, Inv1, A1);   // Inv1 = A^{-1} (mod p)
    if (d ==0) return false;
    conv(tmp, Inv1); // convert to long for translation to a mod-p^r object
  }
  // back to working mod p^r

  Mat<zz_p> XX;
  conv(XX, tmp); // XX = A^{-1} (mod p)

  // Now lift the solution modulo p^r

  // Compute the "correction factor" Z, s.t. XX*A = I - Z (mod p^r)
  long n = A.NumRows();

  const Mat<zz_p> I = ident_mat_zz_p(n);
  Mat<zz_p> T;

  Mat<zz_p> Z; // = I - XX*A
  mul(T, XX, A);
  sub(Z, I, T);

  // The inverse of A is ( I+(Z)+(Z)^2+...+(Z)^{r-1} )*XX (mod p^r). We use
  // O(log r) products to compute it as (I+Z)* (I+(Z)^2)* (I+(Z)^4)*...* XX

  long e = NextPowerOfTwo(r); // 2^e is smallest power of two >= r

  Mat<zz_p> prod; // = I + Z;
  add(prod, I, Z);
  for (long i=1; i<e; i++) {
    sqr(Z, Z);     // = (Z)^{2^i}

    // prod *= (I+Z) = sum_{j=0}^{2^{i+1}-1} (Z)^j
    mul(T, prod, Z);
    add(prod, prod, T);
  }

  mul(X, prod, XX); // X = A^{-1} mod p^r
  return true;
#endif
}




// Preconditions:
//   * MM and M have same dimensions
//   0 <= MM[i][j] < P
//   gcd(p, P) = 1, where p = zz_p::modulus()
//   p*P can be represented as a long
// Postcondition:
//   0 <= MM'[i][j] < P*p
//   MM'[i][j] = MM[i][j] mod P and MM'[i][j] = M[i][j] mod p
//
//  Here, MM' is the new value of MM and MM is the original value

void CRT(Mat<long>& MM, long P, const Mat<zz_p>& M)
{
   long n = MM.NumRows();
   long m = MM.NumCols();

   if (M.NumRows() != n || M.NumCols() != m)
      LogicError("CRT: dimension mismatch");

   if (P == 1) {
      conv(MM, M);
      return;
   }

   long p = zz_p::modulus();
   mulmod_t pinv = zz_p::ModulusInverse();

   sp_reduce_struct p_aux = sp_PrepRem(p);
   // FIXME: NTL should make this visible


   long P_inv = InvMod(rem(cast_unsigned(P), p, p_aux), p);
   mulmod_precon_t P_inv_precon = PrepMulModPrecon(P_inv, p, pinv);

   for (long i = 0; i < n; i++) {
      for (long j = 0; j < m; j++) {
         long x = MM[i][j];
         long y = rep(M[i][j]);


         // compute h = (y-x)*P_inv mod p
         long x_reduced = rem(cast_unsigned(x), p, p_aux);
         long diff = SubMod(y, x_reduced, p);
         long h = MulModPrecon(diff, P_inv, p, P_inv_precon);

         MM[i][j] = x + P*h;
      }
   }
}

// Returns a list of prime factors and their multiplicity,
// N = \prod_i (factors[i].a)^{factors[i].b}
void factorize(Vec< Pair<long,long> > &factors, long N)
{
  factors.SetLength(0);

  if (N < 2) return;

  PrimeSeq s;
  long n = N;
  while (n > 1) {
    if (ProbPrime(n)) { // n itself is a prime, add (n,1) to the list
      append(factors, cons(n, 1L));
      return;
    }

    long p = s.next();
    if ((n % p) == 0) { // p divides n, find its multiplicity
      long e = 1;
      n = n/p;
      while ((n % p) == 0) {
        n = n/p;
        e++;
      }
      append(factors, cons(p, e)); // add (p,e) to the list
    }
  }
}


// VJS: I re-worked this to use zz_p::modulus(), instead of
// passing p as a parameter

// Convert to zz_p and invert mod p, assumes that p is easy to factor
// Returns true on success, false on failure
bool invMod(Mat<zz_p>& X, const Mat<long>& Y)
{
  long p = zz_p::modulus();
  Mat<long> XX(INIT_SIZE, Y.NumRows(), Y.NumCols());
  clear(XX);
  long PP = 1;

  Vec< Pair<long,long> > factors;
  factorize(factors, p);  // prime-power factorization

  { zz_pPush push; // backup modulus

    for (long i=0; i<factors.length(); i++) {
      long pr = power_long(factors[i].a, factors[i].b);
      zz_p::init(pr);
      mat_zz_p Ymod, Xmod;
      conv(Ymod, Y);             // convert to mat_zz_p
      if (!ppInvert(Xmod, Ymod, factors[i].a, factors[i].b)) // invert mod pi^ei
        return false; // inversion failed
      CRT(XX, PP, Xmod);
      PP *= pr;
    }
  }
  // modulus restored

  conv(X, XX); // convert result to mat_zz_p
  return true; // succeeded
}


// Experiment with faster random number generation...maybe we can use this elsewhere,
// like in helib.

#define RBH_BUFSZ (1024)

struct RandomBndHelper {

   long p;
   long nb;
   unsigned long mask;

   RandomStream *str;

   RandomBndHelper() : p(0) { }

   void init(long _p)
   {
      if (_p <= 1) LogicError("RandomBndHelper::init: bad args");

      if (!p) {
         str = &GetCurrentRandomStream();
      }

      p = _p;
      long l = NumBits(p-1);
      nb = (l+7)/8;
      mask = (1UL << l)-1UL;
   }

   long next()
   {
      unsigned char buf[NTL_BITS_PER_LONG/8];
      long tmp;

      do {
         str->get(buf, nb);

         unsigned long word = 0;
         for (long i = nb-1; i >= 0; i--) word = (word << 8) | buf[i];

         tmp = long(word & mask);
      } while (tmp >= p);

      return tmp;
   }
};




void RandomFill(Mat<zz_p>& X)
{
    FHE_TIMER_START;

   long n = X.NumRows();
   long m = X.NumCols();

   RandomBndHelper H;
   H.init(zz_p::modulus());


   for (long i = 0; i < n; i++) {
      zz_p *row = X[i].elts();
      for (long j = 0; j < m; j++) {
         row[j].LoopHole() = H.next();
      }
   }
}

// equivalent to conv(X, A), but assumes entries of A
// are already reduced
void conv_reduced(Mat<zz_p>& X, const Mat<long>& A)
{
   long n = A.NumRows();
   long m = A.NumCols();

   X.SetDims(n, m);
   for (long i = 0; i < n; i++) {
      zz_p *Xrow = X[i].elts();
      const long *Arow = A[i].elts();
      for (long j = 0; j < m; j++)
         Xrow[j].LoopHole() = Arow[j];
   }
}


// generates a random pair of nxn matrices with X*Y = I
// assumes zz_p::modulus() is easy to factor
void GenerateMatrixPair(Mat<zz_p>& X, Mat<zz_p>& Y, long n)
{
  FHE_TIMER_START;

   if (n <= 0) LogicError("GenerateMatrixPair: bad args");

   Mat<long> XX, YY;
   XX.SetDims(n, n);
   clear(XX);

   YY.SetDims(n, n);
   clear(YY);

   long P = 1;

   for (long p = 2, f = zz_p::modulus(); f > 1; p++) {
      if (f % p == 0) {
         long r = 1;
         long pr = p;
         f /= p;

         while (f % p == 0) {
            r++;
            pr *= p;
            f /= p;
         }

         zz_pPush push(pr);

         Mat<zz_p> X1, Y1;
         X1.SetDims(n, n);
         Y1.SetDims(n, n);

         do {
            RandomFill(X1);
         } while (!ppInvert(Y1, X1, p, r));

         CRT(YY, P, Y1);
         CRT(XX, P, X1);
         P *= pr;
      }
   }

   conv_reduced(X, XX);
   conv_reduced(Y, YY);

   FHE_TIMER_STOP;
}

// advance the input stream beyond white spaces and a single instance of cc
void seekPastChar(std::istream& str, int cc)
{
   int c = str.get();
   while (isspace(c)) c = str.get();
   if (c != cc) {
     std::cerr << "Searching for cc='"<<(char)cc<<"' (ascii "<<cc<<")"
	       << ", found c='"<<(char)c<<"' (ascii "<<c<<")\n";
     exit(1);
   }
}


#ifdef TEST_TOOLS
int main()
{
  long n = 200;

  long f = 3L*3L*3L*5L*7L*7L*7L*7L;

  zz_p::init(f);

  Mat<zz_p> X, Y, I;
  GenerateMatrixPair(X, Y, n);
  I = ident_mat_zz_p(n);

  if (X*Y == I)
    cerr << "GOOD1\n";
  else
    cerr << "BAD1\n";

  Mat<long> M(INIT_SIZE, n, n);
  do {
    for (long i=0; i<n; i++) for (long j=0; j<n; j++)
      M[i][j] = RandomBnd(100) - 50; // [-50,50)
  } while (!invMod(Y, M));

  conv(X, M);
  if (X*Y == I)
    cerr << "GOOD2\n";
  else
    cerr << "BAD2\n";
}

#endif // ifdef TEST_TOOLS
