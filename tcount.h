#ifndef TCOUNT_H
#define TCOUNT_H

#include "rint.h"

#include <vector>

#include <iostream>

struct min_unitaries
{
  min_unitaries() : min_t_count(-1), k(0), m(0), factor_calls(0), norm_solver_calls(0) {}

  long min_t_count;
  zwt x;
  std::vector<zwt> y;
  int k; // determinant is \w^k
  long m; // power of \sqrt{2} in denominator
  void to_canonical_form();
  bool operator == ( const min_unitaries& rhs ) const;

  long factor_calls;
  long norm_solver_calls;
};

std::ostream& operator<< (  std::ostream& out , const min_unitaries& mu );
min_unitaries min_t_count( const zwt& x, long m, int k );


#endif // TCOUNT_H
