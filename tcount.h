#ifndef TCOUNT_H
#define TCOUNT_H

#include "rint.h"

#include <vector>

#include <iostream>

struct min_unitaries
{
  long min_t_count;
  zwt x;
  std::vector<zwt> y;
  int k; // determinant is \w^k
  long m; // power of \sqrt{2} in denominator
  void to_canonical_form();
  bool operator == ( const min_unitaries& rhs ) const;
};

std::ostream& operator<< (  std::ostream& out , const min_unitaries& mu );
min_unitaries min_t_count( const zwt& x, long m, int k );


#endif // TCOUNT_H
