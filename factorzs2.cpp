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


#include "factorzs2.h"
#include "appr/normsolver.h"

#include <cassert>

zfactorization factorize(const ztype &val)
{
  return normSolver::instance().factor(val);
}

//factors into zs2 primes
zs2factorization factorize(const zs2type &val, const zfactorization &factors)
{
  zs2factorization res;
  res.solvable = true;

  auto rem = val;
  size_t pos = 0;

  if( factors.prime_factors.size() > 0  &&  factors.prime_factors[0].first == 2 )
  {
    zs2type sol2(2,1);
    for( int i = 0; i < factors.prime_factors[0].second; ++i )
    {
      rem /= sol2;
    }
    res.ramified_prime_power = factors.prime_factors[0].second;
    pos++;
  }

  const auto& ns = normSolver::instance();
  for( ;pos < factors.prime_factors.size(); ++pos )
  {
    const auto& v = factors.prime_factors[pos];
    auto md = v.first % 8;
    if( md == 1 || md == 7 ) // prime splits in Z[\sqrt{2}]
    {
      // solve norm equation in Z[\sqrt{2}] and do proper factorization
      zs2type ans;
      ns.solve(v.first,ans);
      zs2type ans_cnj = ans.g_conjugate();
      int r1 = 0;
      int r2 = 0;
      for( int i = 0; i < v.second; ++i )
      {
        if( rem.divides(ans) )
        {
          rem /= ans;
          r1 ++;
        }
        else
        {
          rem /= ans_cnj;
          r2 ++;
        }
      }

      if( r1 > 0 )
        res.prime_factors.push_back(std::make_pair(ans,r1));

      if( r2 > 0 )
        res.prime_factors.push_back(std::make_pair(ans_cnj,r2));

      if( md == 7 )
      {
        if( (r1 % 2 != 0 ) && (r1 % 2 != 0 ) )
          res.solvable = false;
      }
    }
    else if( md == 3 || md == 5 ) // prime inert in Z[\sqrt{2}]
    {
      for( int i = 0; i < (v.second / 2); ++i )
      {
        rem /= v.first;
      }
      res.prime_factors.push_back(std::make_pair(zs2type(v.first,0),v.second/2));
    }
  }

  assert( rem.norm() == 1 );

  auto ulg = unit_log(rem);
  res.unit_power =ulg.second;
  res.sign = ulg.first;

  return res;
}

// necessary condition only
bool is_solvable(const zfactorization &factors)
{
  //assume that factors.prime_factors is sorted

  if( factors.prime_factors.size() == 0 )
    return factors.sign > 0;

  size_t pos = 0;
  if( factors.prime_factors[0].first == 2 )
    pos++;

  //from this point all primes are odd
  for( ; pos < factors.prime_factors.size(); ++pos )
  {
    const auto& v = factors.prime_factors[pos];
    auto md = v.first % 8;
    if( md == 1 ) // prime splits completely
      continue;
    else
    { // corresponds to the one of quadratic subrings Z[i], Z[\sqrt{2}], Z[\sqrt{-2}]
      if( v.second % 2 == 0 )
        continue;
      else
        return false;
    }
  }
  return true;
}
