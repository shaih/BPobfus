#ifndef VEC_L_H_INCLUDED
#define VEC_L_H_INCLUDED
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
/* vec_l - functions for  a vector of type long */

#include <NTL/tools.h>
#include <NTL/vector.h>

typedef NTL::Vec<long> vec_l;

void sub(vec_l& x, const vec_l& a, const vec_l& b);
void mul(long& res, const vec_l& x, const vec_l& a);

inline vec_l& operator-=(vec_l& x, const vec_l& a)
{ sub(x, x, a); return x; }

vec_l operator-(const vec_l& a, const long& b);
vec_l operator-(const vec_l& a, const vec_l& b);

#endif // VEC_L_H_INCLUDED
