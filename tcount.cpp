#include "tcount.h"
#include "factorzs2.h"
#include "solvenormequation.h"
#include "matrix2x2.h"

#include "output.h"

#include <iostream>
#include <cassert>

using namespace std;

std::pair<int,long> min_w_t_count( const zwt& x, const zwt& y, long m, long k )
{
  zwt ycm = -y.conjugate();
  matrix2x2<mpz_class> mtr(x, ycm * zwt::omega(k),
                           y, x.conjugate() * zwt::omega(k), m );
  assert(mtr.is_unitary());

  int min_i = 0;
  long min_tc = mtr.t();

  for( int i = 1; i < 8; ++i )
  {
    mtr.d[1][0] = y * zwt::omega(i);
    mtr.d[0][1] = ycm * zwt::omega(8-i+k);
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
  res.x = x;
  res.m = m;
  res.k = k;

  std::vector<zwt> y[2];

  long n0 = 2*m - x.abs2().gde();
  assert( n0 >= 4 );

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
  }

  if( y[0].empty() )
  {
    res.min_t_count = n0 - 1;
    swap( y[1], res.y );
  }
  else
  {
    res.min_t_count = n0 - 2;
    swap( y[0], res.y );
  }


  return res;
}
