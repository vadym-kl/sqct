#ifndef SOLVENORMEQUATION_H
#define SOLVENORMEQUATION_H

#include "factorzs2.h"

struct norm_equation_solution
{
  bool exists;
  std::vector< std::pair<zwt,long> > ramified;
  std::vector< std::pair<zwt,long> > split;
  std::vector< std::pair<zwt,long> > inert;
};

norm_equation_solution solve_norm_equation( const zs2type& rhs );
std::vector<zwt> all_solutions( const norm_equation_solution& neq_sln );

#endif // SOLVENORMEQUATION_H
