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


#include "tcount.h"
#include "factorzs2.h"
#include "solvenormequation.h"
#include "matrix2x2.h"

#include "output.h"

#include <iostream>
#include <algorithm>
#include <cassert>

using namespace std;

template < class T >
int chop( const T& x )
{
  if( x > 0 )
  {
    return mpz_class(x & 0xFFFF).get_si();
  }
  else
    return -mpz_class((-x) & 0xFFFF).get_si();
}

template < class T >
ring_int<int> chop( const ring_int<T>& x )
{
  ring_int<int> r(chop(x.v[0]),chop(x.v[1]),chop(x.v[2]),chop(x.v[3]));
  return r;
}

std::pair<int,long> min_w_t_count( const zwt& xc, const zwt& yc, long m, long k )
{
  ring_int<int> x = chop(xc);
  ring_int<int> y = chop(yc);

  ring_int<int> ycm = -y.conjugate();
  matrix2x2<int> mtr(x, ycm * ring_int<int>::omega(k),
                           y, x.conjugate() * ring_int<int>::omega(k), m );
  //assert(mtr.is_unitary());

  int min_i = 0;
  long min_tc = mtr.t();

  for( int i = 1; i < 2; ++i )
  {
    mtr.d[1][0] = y * ring_int<int>::omega(i);
    mtr.d[0][1] = ycm * ring_int<int>::omega(8-i+k);
    int t = mtr.t();
    if( t < min_tc )
    {
      min_tc = t;
      min_i = i;
    }
  }
  return std::make_pair(min_i,min_tc);
}

min_unitaries min_t_count(const zwt &x, long m, int k)
{
  min_unitaries res;
  res.min_t_count = -1;
  res.x = x;
  res.m = m;
  res.k = k;

  std::vector<zwt> y[4];

  long n0 = 2*m - x.abs2().gde();
//  if( n0 < 4 )
//  {
//    cerr << x << "," << m << "," << k << endl;
//    throw std::logic_error(":unsupported regime");
//  }

  zs2type pow2 = zs2type(ztype(1) << m,0);
  zs2type rhs(pow2 - x.abs2());

  auto sln = solve_norm_equation(rhs);
  if( sln.exists )
  {
    auto all_slns = all_solutions(sln);
    for( const auto& a : all_slns )
    {
      auto tc = min_w_t_count(x,a,m,k);
      int d = tc.second - n0 + 2;
      assert( d == 0 || d == 1 );
      y[d].push_back(a * zwt::omega(tc.first));
    }

    for( int d = 0; d < 4 ;++d )
    {
      if( ! y[d].empty() )
      {
        res.min_t_count = n0 - 2 + d;
        swap( y[d], res.y );
        //res.to_canonical_form();
        res.factor_calls = sln.factor_calls;
        res.norm_solver_calls = sln.norm_solver_calls;
        return res;
      }
    }

  }

  res.factor_calls = sln.factor_calls;
  res.norm_solver_calls = sln.norm_solver_calls;
  return res;
}


void min_unitaries::to_canonical_form()
{
  if( x < -x )
    x = -x;

  for( size_t i = 0; i < y.size(); ++i )
  {
    y[i] = y[i].i_canonical();
  }

  {
    vector<zwt> tmp;
    tmp.resize(y.size());
    sort(y.begin(),y.end());
    auto r = unique_copy(y.begin(),y.end(),tmp.begin());
    tmp.resize(distance(tmp.begin(),r));
    swap(tmp,y);
  }

  for( size_t i = 0; i < y.size(); ++i )
  {

    matrix2x2<mpz_class> mm1(x   ,-y[i].conjugate() * zwt::omega(k),
                      y[i], x.conjugate() * zwt::omega(k)    , m);

    matrix2x2<mpz_class> mm2(x                    ,-y[i].conjugate() * zwt::omega(k+7),
                       y[i] * zwt::omega(1) , x.conjugate() * zwt::omega(k)     , m);

    if( mm1.t() == mm2.t() )
    {
      y[i] = y[i].w_canonical();
    }
  }

  {
    vector<zwt> tmp;
    tmp.resize(y.size());
    sort(y.begin(),y.end());
    auto r = unique_copy(y.begin(),y.end(),tmp.begin());
    tmp.resize(distance(tmp.begin(),r));
    swap(tmp,y);
  }

}

bool min_unitaries::operator ==(const min_unitaries &rhs) const
{
  if( y.size() == 0 && rhs.y.size() == 0 ) //empty results are equal
    return true;

  if( rhs.k != k ) return false;
  if( rhs.m != m ) return false;
  if( rhs.min_t_count != min_t_count ) return false;
  if( rhs.x != x ) return false;
  if( rhs.y.size() != y.size() ) return false;

  for( size_t s = 0; s < y.size(); ++s )
  {
    if( y[s] != rhs.y[s] )
      return false;
  }

  return true;
}

string min_unitaries::short_title()
{
  stringstream ss;
  ss << "{\"x[0]\",\"x[1]\",\"x[2]\",\"x[3]\"}" << "," <<
        "{\"y[0]\",\"y[1]\",\"y[2]\",\"y[3]\"}" << "," <<
        "\"m\"," <<
        "\"k\"," << "\"Solutions\"";
  return ss.str();
}

string min_unitaries::short_str() const
{
  stringstream ss;
  if( y.size() > 0 )
    ss << x << "," << y[0] << "," << m << "," << k << "," << y.size();
  else
    ss << x << "," <<  "{0,0,0,0}" << "," << m << "," << k << "," << y.size();
  return ss.str();
}

min_unitaries::operator matrix2x2<mpz_class>() const
{
  matrix2x2<mpz_class> r(x,-y[0].conjugate() * zwt::omega(k),
      y[0] , x.conjugate()*zwt::omega(k), m);
  return r;
}


ostream &operator <<(ostream &out, const min_unitaries &mu)
{
  out << "{" << mu.k << "," << mu.m <<
          "," << mu.x <<
          "," <<  mu.min_t_count;

  out << ",{";
  if(mu.y.size() > 0 )
    cout << mu.y[0];

  for( size_t i = 1; i < mu.y.size(); ++i )
  {
    out << "," << mu.y[i];
  }
  out << "}}";
  return out;
}
