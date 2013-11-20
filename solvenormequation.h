#ifndef SOLVENORMEQUATION_H
#define SOLVENORMEQUATION_H

#include "factorzs2.h"

struct norm_equation_solution
{
  norm_equation_solution() : exists(false), unit_power(0), norm_solver_calls(0), factor_calls(0) {}
  bool exists;
  long unit_power;
  std::vector< std::pair<zwt,long> > ramified;
  std::vector< std::pair<zwt,long> > split;
  std::vector< std::pair<zwt,long> > inert;

  long norm_solver_calls;
  long factor_calls;

  operator zwt()
  {
    zwt r(1,0,0,0);
    for( const auto& a : ramified )
    {
      for( long i = 0; i < a.second; ++i )
      {
        r = r * a.first;
      }
    }

    for( const auto& a : split )
    {
      for( long i = 0; i < a.second; ++i )
      {
        r = r * a.first;
      }
    }

    for( const auto& a : inert )
    {
      for( long i = 0; i < a.second; ++i )
      {
        r = r * a.first;
      }
    }

    r = r * ::unit_power<mpz_class>(std::make_pair(1,unit_power));

    return r;
  }
};

norm_equation_solution solve_norm_equation( const zs2type& rhs );
std::vector<zwt> all_solutions( const norm_equation_solution& neq_sln );

#endif // SOLVENORMEQUATION_H
