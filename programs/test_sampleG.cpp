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
/**************************************************************
 * A program for testing the sampleG module
 **************************************************************/
//#include <mpfr.h>
#include <NTL/lzz_p.h>

NTL_CLIENT
#include "../utils/argmap.h"
#include "../TDMatrix.h"

#define NFACTORS 8
long tinyFactors[NFACTORS] = { 2, 3, 5, 7, 11, 13, 17, 19 };

int main(int argc, char *argv[])
{
#ifdef DEBUG
  NTL::SetSeed((unsigned char*)"test_sampleG",12);
#endif
  ArgMapping amap;

  long k = 2;
  amap.arg("k", k, "number of factors (~ log q)");

  long n = 2;
  amap.arg("n", n, "the small dimension"); // no default info

  long e=2;
  amap.arg("e", e, "the power of each factor");

  amap.parse(argc, argv); // parses and overrides initail values
                         // for each parameter, one per line,

  assert (k <= NFACTORS); // for this test we don't want many factors
  long m = (2+ k*e)*n;

  TDMatrixParams params(n,k,e,m);
  Vec<long>& factors = params.factors;

  cout << "n="<<n<<", m="<<m<<", k="<<k<<", e="<<e<<", q="<<params.getQ()<<endl;

  Vec<vec_zz_p> syndrome(INIT_SIZE, k);
  for (long i=0; i<k; i++) {
    params.zzp_context[i].restore();
    syndrome[i].SetLength(n);
    for (long j=0; j<n; j++) syndrome[i][j] = random_zz_p();
    cout << "syndrome mod "
	 << power_long(factors[i],e)<<"= "<<syndrome[i]<<endl <<flush;
  }

  // Sample x such that G*x = syndrome
  Vec<long> x;
  sampleG(x, syndrome, &params);

  cout << "sampled vector="<<x<<endl;
  assert(x.length()== n*e*k);
  // check that we have the right answer modulo p^e for all factors

  for (long i=0; i<k; i++) {
    params.zzp_context[i].restore();
    vec_zz_p xMod = conv<vec_zz_p>(x);
    long index = 0;
    for (long j=0; j<n; j++) { // check the next syndrome entry mod f_i
      zz_p xp = to_zz_p(0L);
      zz_p curFactor = to_zz_p(1L);
      for (long f=0; f<k; f++) for (long ei=0; ei<e; ei++) {
	  xp += (curFactor * x[index++]);
	  curFactor *= factors[f];
      }
      assert(syndrome[i][j] == xp);
    }
  }
  cout << "PASSED\n";

  return 0;
}
