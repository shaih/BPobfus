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
/*********************************************************************
vec_l: A module for handling vectors of type long
**********************************************************************/
#include <NTL/matrix.h>
#include "mat_l.h"
#include "vec_l.h"

#include <NTL/new.h>
#include <NTL/vec_long.h>

NTL_CLIENT

//x = a-b
void sub(vec_l& x, const vec_l& a, const vec_l& b)
{
    long n = a.length();
    if (b.length() != n) LogicError("vector sub: dimension mismatch");
    x.SetLength(n);

    for (long i = 0; i < n; i++)
        x[i] = a[i] - b[i];

}

//res = x*a, naive implementation

void mul(long& res, const vec_l& x, const vec_l& a)
{
    long val;
    long n = a.length();
    long i;
    val = 0;
    for (i = 0; i < n; i++)
    {
        val += a.at(i)*x.at(i);
    }
    res = val;
}

//res = a-b
vec_l operator-(const vec_l& a, const vec_l& b)
{
    vec_l res;
    sub(res, a, b);
    return res;
    //NTL_OPT_RETURN(vec_l, res);
}

//res = a-b
vec_l operator-(const vec_l& a, const long& b)
{
    vec_l res;

    int n = a.length();
    int i;
    res.SetLength(n);
    for (i = 0; i< n; i++)
        //res.put(i,(a.get(i) - b));
        res[i] = a[i] - b;
    return res;
    //NTL_OPT_RETURN(vec_l, res);
}

