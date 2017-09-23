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
/////////////////////////////
// sampleG - Sampling pre-image z such that G*z = v, G is the "easy" matrix
//////////////////////////////////////////

//#include <mpfr.h>
#include <NTL/mat_lzz_p.h>
#include "utils/timing.h"

NTL_CLIENT

#include "mat_l.h"
#include "stash.h"
#include "TDMatrix.h"

static zz_p
updateSyndromValue(vec_l &xOut, const Vec<vec_zz_p>& syndrom,
		   TDMatrixParams* params, long lFactor, long row,long offset);

//#define DEBUG

// Sample from a discrete Gaussian until you get a sample point equal
// to u modulo f. Uses a stack to put points that are not equal to u.
static long sampleMod(long u, long f,
                      Gaussian1Dsampler& sampler, SampleStash& stash)
{
  FHE_TIMER_START;

  long umodf = u%f;
  if (umodf<0) umodf += f;

  long val = (long)stash.getPoint(umodf);
  if (val != umodf+1) // found in stash
    return val;

  // If not found in stash, sample until you get it
  while (true) {
    NTL::Pair<bool,long> extra(false,0);
    val = sampler.getSample(/*mu=*/0.0, extra);
    long valmodf = val%f;
    if (valmodf<0) valmodf += f;

    if ( valmodf == umodf ) // found
      return val;
    else                    // not found, keep point for later
      stash.setPoint(valmodf, (int)val);
  }
}


// SampleG: Sample vector xOut such that G * xOut = syndrom. syndrom is given
// in CRT representation but x is a vector of longs.

// Generate a vector of size params>n*kFactors*e
// For each row in G, choose kFactors*e elements of xOut

// The function updates the syndrom vector, based on all previously
// chosen values in xOut. We update syndrom[lFactor][row]
// only just before we need to use it (as opposed to updating
// it for every new value of xOut as it is chosen)

// based on the updated syndrome, we sample a new entry in xOut,
// equal to uVal mod factor and append it to xOut

int sampleG(vec_l &xOut,
            const Vec<vec_zz_p>& syndrom, TDMatrixParams* params)
{
  FHE_TIMER_START;
  xOut.SetLength(0);
  xOut.SetMaxLength(params->m); // allocate space

  long sigma = params->r * params->maxFactor;
  Gaussian1Dsampler sampler(sigma);

  zz_pPush push; // backup the NTL current modulus

  // For each row in G, choose kFactors*e elements of xOut
  for (long row = 0; row < params->n; row++)
    {
      long offset = row * (params->kFactors) * (params->e);

      for (long lFactor = 0; lFactor < params->kFactors; lFactor++)
        {
	  long factor = params->factors[lFactor];
	  params->zzp_context[lFactor].restore(); // NTL-modulus := factor

	  // lazy update of the syndrom vector, based on all previously
	  // chosen values in xOut. We update syndrom[lFactor][row]
	  // only just before we need to use it (as opposed to updating
	  // it for every new value of xOut as it is chosen).

	  zz_p zzNewVal = updateSyndromValue(xOut, syndrom,params, lFactor, row, offset);

	  long uVal = conv<long>(zzNewVal);

	  // Done updating the uSyndrom, now choose next e elements in xOut

	  for (long lPower = 0; lPower < params->e; lPower++)
            {
	      // sample a new entry in xOut, equal to uVal mod factor
	      long samp = sampleMod(uVal, factor, sampler, params->stash[lFactor]);
	      assert( 0 == ((uVal-samp)%factor) ); // sanity check

	      xOut.append(samp);

	      // update uVal := (uVal-samp)/factor (over the integers)
	      uVal -= samp;
	      uVal /= factor;
            }
        }
    }
  return 0;
}


// lazy update of the syndrom vector, based on all previously chosen values
// in xOut. We update syndrom[lFactor][row] just before we need to use it
// (as opposed to updating it for every new value of xOut as it is chosen).

static zz_p
updateSyndromValue(vec_l &xOut, const Vec<vec_zz_p>& syndrom,
		   TDMatrixParams* params, long lFactor, long row, long offset)
{
    zz_p zzNewVal = syndrom[lFactor][row];
#ifdef DEBUG
    cout << zzNewVal << endl;
#endif
    long xIndex=offset;
    for (long vFactor = 0; vFactor < lFactor; vFactor++)
        for (long ie=0; ie < params->e; ie++) // repeat e times
        {
            long xVal = xOut[xIndex++];     // xIndex = offset+ie
            zzNewVal -= xVal;
            zzNewVal *= params->fInv[vFactor][lFactor];
            // update via u := (u - xval)/pi (mod pj)
#ifdef DEBUG
            cout << zzNewVal << endl;
#endif
        }
    return zzNewVal;
}
