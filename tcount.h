#ifndef TCOUNT_H
#define TCOUNT_H

#include "rint.h"

#include <vector>



struct min_unitaries
{
  long min_t_count;
  zwt x;
  std::vector<zwt> y;
  int k; // determinant is \w^k
  long m; // power of \sqrt{2} in denominator
};

min_unitaries min_t_count( const zwt& x, long m, int k );


#endif // TCOUNT_H
