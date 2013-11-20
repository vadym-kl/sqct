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


#include "rcup.h"
#include "findhalves.h"

#include <array>
#include <algorithm>

#include "output.h"

#include "timemeasurement.h"

using namespace std;

double min_positive_first( const halves_t& hv )
{
  halves_t::value_type z(0,0);
  auto r = std::upper_bound(hv.begin(),hv.end(),z);
  if( r != hv.end() )
    return r->first;
  else
    return 1.;
}

zwt to_zwt( const std::array<long,5>& val )
{
  return zwt(val[0],val[1]+val[3],val[2],val[3]-val[1]);
}

typedef std::pair< halves_t::const_iterator, halves_t::const_iterator> it_interval_t;

it_interval_t get_interval(const halves_t& hv, double min_v, double max_v )
{
  it_interval_t r;
  halves_t::value_type l(min_v,0);
  halves_t::value_type h(max_v,0);
  r.first = std::upper_bound(hv.begin(),hv.end(),l);
  if( r.first != hv.begin() ) --r.first;
  r.second = std::upper_bound(hv.begin(),hv.end(),h);
  return r;
}

void rcup::merge_halves( const halves_t& re, const halves_t& im, int k )
{
  size_t sz = re.size();

  for( size_t j = 0; j < sz; ++j )
  {
    const auto& re_j = re[j];

    auto II = get_interval( im, I.first - re_j.first,I.second - re_j.first );
    for( auto jj = II.first; jj != II.second; ++jj )
    {
      double val = jj->first + re_j.first;
      if( val >= IL.first && val <= IL.second )
      {
        long c = recover_imag(*jj,k);
        long a = recover_real(re_j,k);
        out.push_back(make_pair(jj->first + re_j.first,std::array<long,5>{a,re_j.second,c,jj->second,k}));
      }
    }
  }
}

long rcup::recover_real(const halves_t::value_type &a, int k) const
{
  hprr v = cos2m[k] - a.second * hprHelpers::sqrt2();
  hprr v_minus_a = hprr(a.first) / reW[k];
  return to_long(v - v_minus_a);
}

long rcup::recover_imag(const halves_t::value_type &a, int k) const
{
  hprr v = sin2m[k] - a.second * hprHelpers::sqrt2();
  hprr v_minus_a = hprr(a.first) / imW[k];
  return to_long(v - v_minus_a);
}

rcup::rcup(long n, const hprr &phi, const hprr &delta)
{
  obs.n = n; obs.phi = phi; obs.delta = delta;

  R.first = delta;

  long m = floor(n/2) + 3;
  halves_t L_re[2];
  halves_t L_im[2];
  for( int k = 0; k < 2; ++k )
  {
    hprr theta = - (phi / hprr(2) -hprHelpers::pi()*k/hprr(8));
    cos2m[k] = cos( theta ) * sqrt2pow(m);
    sin2m[k] = sin( theta ) * sqrt2pow(m);
    reW[k] = weight( cos( theta ) ,m );
    imW[k] = weight( sin( theta ) ,m );

    timepoint t1,t2;
    timepoint t1r,t2r;

    get_timepoint(t1);
    get_timepoint_real(t1r);
    L_re[k] = findhalves3(cos( theta ), m, delta  );
    L_im[k] = findhalves3(sin( theta ), m, delta  );
    get_timepoint(t2);
    get_timepoint_real(t2r);

    obs.halves_size += (L_re[k].size() + L_im[k].size());
    obs.find_halves_real_time += time_diff_d(t1r,t2r);
    obs.find_halves_cpu_time += time_diff_d(t1,t2);
  }

  I = interval_t( 0,min( min_positive_first(L_re[0]), min_positive_first(L_re[1]) ) );
  // Question 1: is doulbe precision sufficient for the purpose of this part ? -- assume yes
  // Question 2: how likely is to get repetitions of epsilon -- assume that it is not likely

  double dl = to_ld(delta)*to_ld(delta);
  I.second = min(dl,I.second);

  while( I.first < dl )
  {
    IL.first = I.first * (1. - 1e-5 );
    IL.second = I.second * (1. + 1e-5 );

    timepoint t1,t2;
    get_timepoint(t1);

    for( int k = 0; k < 2; ++k )
      merge_halves( L_re[k], L_im[k], k);

    get_timepoint(t2);
    obs.merge_halves_time += time_diff_d(t1,t2);

    std::sort(out.begin(),out.end());

    obs.max_tuples_memory_size = std::max( obs.max_tuples_memory_size, out.size() );
    obs.tuples_processed += out.size();

    size_t sz_max = out.size();

    timepoint t3,t4;
    get_timepoint(t3);

    for( size_t j = 0; j < sz_max; ++j )
    {
      const auto& o = out[j];
      zwt xp = to_zwt(o.second);
      int dm = xp.gde();
      xp.div_eq_sqrt2(dm);

      R.second = min_t_count(xp,m-dm,o.second[4]);
      obs.factor_calls += R.second.factor_calls;
      obs.norm_equation_calls += R.second.norm_solver_calls;
      if( n == R.second.min_t_count )
      {
//        cout << R.second.min_t_count << endl;
//        cout << o.first << endl;
//        cout << xp << endl;
//        cout << m-dm << endl;
//        cout << o.second[4] << endl;
        R.first = sqrt(o.first);
        get_timepoint(t4);
        obs.tcount_time += time_diff_d(t3,t4);
        return;
      }
    }

    get_timepoint(t4);
    obs.tcount_time += time_diff_d(t3,t4);

    out.clear();
    double t = (I.second - I.first) * 2.0;
    I.first = I.second;
    I.second = min(I.first + t,dl);
  }

  R.second.k = 0;
  R.second.m = 0;
  R.second.min_t_count = -1;
  R.second.x = zwt(0,0,0,0);
  R.second.y.clear();
}

