//     Copyright (c) 2013 Vadym Kliuchnikov sqct(dot)software(at)gmail(dot)com
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
//     Based on "Practical approximation of single-qubit unitaries by single-qubit quantum Clifford and T circuits"
//     by Vadym Kliuchnikov, Dmitri Maslov, Michele Mosca
//     http://arxiv.org/abs/1212.6964v1, new version of the paper describing this modification: [to be published]

#include "approxlist.h"

#include <stdexcept>
#include <gmpxx.h>

using namespace std;

typedef __int128 ul;

static long calc_b_min( const approxlist_params& in)
{
  hprr tmp = max( -pow2(in.n),
             ( pow2(in.n) * (in.s - hprr(1.)) - hprr(0.5) ) / sqrt( hprr(2.) ) );

  if( !mpfr_fits_slong_p(tmp._x,MPFR_RNDD) )
    throw logic_error("too big search space");

  return mpfr_get_si(tmp._x,MPFR_RNDD);
}

static long calc_b_max( const approxlist_params& in )
{
  hprr tmp = min( pow2(in.n),
               ( pow2(in.n) * (in.s + hprr(1.)) + hprr(0.5) ) / sqrt( hprr(2.) ) );

  if( !mpfr_fits_slong_p(tmp._x,MPFR_RNDU) )
    throw logic_error("too big search space");

  return mpfr_get_si(tmp._x ,MPFR_RNDU);
}

bool approx_list_builder::norm_condition( long a, long b )
{
  ul aa = a;
  aa *= a;
  ul b2b = b;
  b2b *= b;
  b2b <<= 1;
  return (b2b + aa) <= pow4n;
}

approx_list_builder::approx_list_builder(const approxlist_params &_in, approx_list &_out) :
  in(_in),out(_out)
{
  if( in.delta > 0.5 || in.delta < 0.0 )
  {
    cout << "delta:" << in.delta << endl;
    throw logic_error("Unsupported value of delta");
  }

  pow4n = 1; pow4n <<= (in.n * 2);

  out.clear();
  long b_min = calc_b_min(in);
  long b_max = calc_b_max(in);

  hprr sqrt2     ( sqrt(hprr(2)) );
  hprr sqrt2st   ( sqrt(hprr(2)) * hprr(in.step) );
  hprr s2n       ( in.s * pow2(in.n) );
  hprr t         ( s2n - hprr(b_min + in.init) * sqrt2 );
  long double sld = to_ld(in.s);

  hprr r = 0.;

  int st = in.step;
  for( long b = (b_min + in.init) ; b <= b_max ; b+=st )
  {
    long a = mpfr_get_si(t._x, MPFR_RNDN);
    r = t - a;
    if( abs(r) < in.delta && norm_condition(a,b) )
      out.push_back( make_pair(to_ld(r) *sld,b) );
    t -= sqrt2st;
  }
}
