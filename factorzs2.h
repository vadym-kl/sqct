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
  std::vector< std::pair<zs2type, long> > prime_factors;
  int sign;
  int unit_power;
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
