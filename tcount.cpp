#include "tcount.h"
#include "factorzs2.h"
#include "solvenormequation.h"

min_unitaries min_t_count(const zwt &x, long m, int k)
{
  min_unitaries res;
  res.x = x;
  res.m = m;
  res.k = k;

  zs2type pow2 = zs2type(ztype(1) << m,0);
  zs2type rhs(pow2 - x.abs2());
  auto sln = solve_norm_equation(rhs);
  if( sln.exists )
  {
    auto all_slns = all_solutions(sln);

  }
  return res;
}
