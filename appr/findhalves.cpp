//     Copyright (c) 2012 Vadym Kliuchnikov sqct(dot)software(at)gmail(dot)com
//
//     This file is part of SQCT.
//
//     SQCT is free software: you can redistribute it and/or modify
//     it under the terms of the GNU Lesser General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.
//
//     SQCT is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU Lesser General Public License for more details.
//
//     You should have received a copy of the GNU Lesser General Public License
//     along with SQCT.  If not, see <http://www.gnu.org/licenses/>.
//

#include "findhalves.h"
#include "fixedpoint.h"
#include <algorithm>
#include <parallel/algorithm>

#include <omp.h>

typedef __int128 ul;

inline bool norm_condition( long a, long b, ul pow2m )
{
  ul aa = a;
  aa *= a;
  ul b2b = b;
  b2b *= b;
  b2b <<= 1;
  return (b2b + aa) <= pow2m;
}

typedef hprr ft;

using namespace std;

hprr weight( const hprr& alpha, int m )
{
  hprr s2m = sqrt2pow(m);
  return alpha / s2m;
}

halves_t findhalves(const hprr& alpha, int m, const hprr& delta)
{
//  cout << fixed;
//  cout.precision(100);
//  cout << alpha << "," << m << "," << delta << endl;

  ul pow2m = (ul(1) << m);

  hprr s2m = sqrt2pow(m);
  hprr W(weight(alpha,m));
  hprr epsilon = delta * s2m * hprHelpers::sqrt2();
  long b = to_long(floor(-s2m));
  long b_max = to_long(ceil(s2m));
  hprr V = alpha * s2m - hprr(b) * hprHelpers::sqrt2();

  halves_t R;

  ft eps(epsilon);
  long double eld(to_ld(epsilon));
  double w = to_double(W);
  ft v(V);
  ft ms2(-hprHelpers::sqrt2());

  if( to_ld(epsilon) < 0.5 ) // optimized version
  {
    while( b <= b_max )
    {
      long a = to_long(round(v));
      long double dst = to_double(v-a);
      if( abs(dst) < eld && norm_condition(a,b,pow2m) )
        R.push_back(make_pair(w*dst,b));
      //cout << alpha << "," << a << "," << b << "," << v << "," << w*to_double(v-a) << endl;
      b++;
      v+=ms2;


    }
  }
  else //generic version
  {
    while( b <= b_max )
    {
      long a_min = to_long(ceil(v-eps));
      long a_max = to_long(floor(v+eps));
      for( long a = a_min; a <= a_max; ++a )
        R.push_back(make_pair(w*to_double(v-a),b));
      b++;
      v+=ms2;
    }

  }

  std::sort(R.begin(),R.end());
  return R;
}


//////////// optimized version //////////////////////////////////

halves_t findhalves2(const hprr& alpha, int m, const hprr& delta)
{
  hprr s2m = sqrt2pow(m);
  hprr W(weight(alpha,m));
  hprr epsilon = delta * s2m * hprHelpers::sqrt2();
  long b = to_long(floor(-s2m));
  long b_max = to_long(ceil(s2m));
  hprr V = alpha * s2m - hprr(b) * hprHelpers::sqrt2();
  ul pow2m = (ul(1) << m);
  hprr ms2 = -hprHelpers::sqrt2();

  halves_t R;

  long double w = to_ld(W);

  if( to_ld(epsilon) < 0.5 ) // optimized version
  {
    grid_iterator gi(V,hprHelpers::sqrt2(),epsilon,m);
    while( b <= b_max )
    {
      if( gi.down_close() )
      {
        if( norm_condition(gi.down_a(),b,pow2m) )
        {
          R.push_back(make_pair(w*gi.down_frac_part(),b));
        }
      }
      else if( gi.up_close() )
      {
        if( norm_condition(gi.up_a(),b,pow2m) )
        {
          R.push_back(make_pair(w*gi.up_frac_part(),b));
        }
      }

      ++b;
      ++gi;
    }
  }
  else //generic version
  {
    while( b <= b_max )
    {
      long a_min = to_long(ceil(V-epsilon));
      long a_max = to_long(floor(V+epsilon));
      for( long a = a_min; a <= a_max; ++a )
        R.push_back(make_pair(w*to_double(V-a),b));
      b++;
      V+=ms2;
    }

  }

  std::sort(R.begin(),R.end());
  return R;
}

//////////// parallel version //////////////////////////////////
halves_t findhalves3_basic(const hprr& alpha, int m, const hprr& delta, int step, int init);

halves_t findhalves3(const hprr& alpha, int m, const hprr& delta)
{
//  cout << fixed;
//  cout.precision(100);
//  cout << alpha << "," << m << "," << delta << endl;

  int threads = omp_get_max_threads();
  std::vector< halves_t > res(threads);
//  cout << "Using " << threads << " threads" << endl;

  #pragma omp parallel
  {
    int thid = omp_get_thread_num();
    res[thid] = findhalves3_basic(alpha,m,delta,threads,thid);
  }

  std::vector< size_t > szs(threads);
  szs[0] = 0;
  for( int i = 1; i < threads; ++i )
    szs[i] = szs[i-1] + res[i-1].size();

  halves_t r(szs[threads-1]+res[threads-1].size());

  #pragma omp parallel
  {
    int thid = omp_get_thread_num();
    std::copy(res[thid].begin(),res[thid].end(),r.begin()+szs[thid]);
  }

  __gnu_parallel::sort(r.begin(),r.end());

  return r;
}

halves_t findhalves3_basic(const hprr& alpha, int m, const hprr& delta, int step, int init)
{
  hprr s2m = sqrt2pow(m);
  hprr W(weight(alpha,m));
  hprr epsilon = delta * s2m * hprHelpers::sqrt2();
  long b = to_long(floor(-s2m) + init);
  long b_max = to_long(ceil(s2m));
  hprr V = alpha * s2m - hprr(b) * hprHelpers::sqrt2();
  ul pow2m = (ul(1) << m);
  hprr ms2 = -hprHelpers::sqrt2() * hprr(step);

  halves_t R;

  long double w = to_ld(W);

  if( to_ld(epsilon) < 0.5 ) // optimized version
  {
    grid_iterator gi(V,hprHelpers::sqrt2()*step,epsilon,m);
    while( b <= b_max )
    {
      if( gi.down_close() )
      {
        if( norm_condition(gi.down_a(),b,pow2m) )
        {
          R.push_back(make_pair(w*gi.down_frac_part(),b));
        }
      }
      else if( gi.up_close() )
      {
        if( norm_condition(gi.up_a(),b,pow2m) )
        {
          R.push_back(make_pair(w*gi.up_frac_part(),b));
        }
      }

      b+=step;
      ++gi;
    }
  }
  else //generic version
  {
    while( b <= b_max )
    {
      long a_min = to_long(ceil(V-epsilon));
      long a_max = to_long(floor(V+epsilon));
      for( long a = a_min; a <= a_max; ++a )
        R.push_back(make_pair(w*to_double(V-a),b));
      b+=step;
      V+= (ms2);
    }

  }

  //std::sort(R.begin(),R.end());
  return R;
}
