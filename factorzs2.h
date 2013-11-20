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


#ifndef FACORZS2_H
#define FACORZS2_H

#include "rint.h"

#include <vector>



struct zfactorization
{
  std::vector< std::pair<ztype, long> > prime_factors;
  int sign;
};

struct zs2factorization
{
  zs2factorization() : solvable(false), ramified_prime_power(0), sign(1), unit_power(0) {}
  bool solvable;
  int ramified_prime_power ;
  std::vector< std::pair<zs2type, long> > prime_factors;
  int sign;
  long unit_power;

  operator zs2type()
  {
    zs2type r(::unit_power<mpz_class>(std::make_pair(sign,unit_power)));
    for( const auto& a : prime_factors )
    {
      for( int k = 0; k < a.second; ++k )
        r = r * a.first;
    }

    zs2type rm(2,1);
    for( int i = 0; i < ramified_prime_power; ++i )
      r = r * rm;

    return r;
  }
};

zfactorization factorize( const ztype& val );
zs2factorization factorize( const zs2type& val, const zfactorization& factors );

bool is_solvable( const zfactorization& factors );

static zs2factorization factorize( const zs2type& val )
{
  zfactorization zf = factorize(val[0]*val[0]-val[1]*val[1]*2);
  return factorize(val,zf);
}


#endif // FACORZS2_H
