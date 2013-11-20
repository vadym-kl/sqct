#include "observations.h"

#include <sstream>

using namespace std;

std::ostream &operator <<(std::ostream &out, const observations &o)
{
  out << std::fixed;
  out.precision(20);

  out << o.n << "," <<
      to_ld(o.phi) << "," <<
      to_ld(o.delta) << "," <<

      o.find_halves_time << "," <<
      o.merge_halves_time << "," <<
      o.tcount_time << "," <<

      o.halves_size << "," <<
      o.b_max << "," <<

      o.max_tuples_memory_size << "," <<
      o.tuples_processed << "," <<
      o.factor_calls << "," <<
      o.norm_equation_calls;
}


observations::observations() :
  n(0),
  phi(0),
  delta(0),

  find_halves_time(0),
  merge_halves_time(0),
  tcount_time(0),

  halves_size(0),
  b_max(0),

  max_tuples_memory_size(0),
  tuples_processed(0),
  factor_calls(0),
  norm_equation_calls(0)
{
}

std::string observations::title()
{
  stringstream out;

  out << "\"n\"" << "," <<
      "\"phi\"" << "," <<
      "\"delta\"" << "," <<

      "\"find_halves_time\"" << "," <<
      "\"merge_halves_time\"" << "," <<
      "\"tcount_time\"" << "," <<

      "\"halves_size\"" << "," <<
      "\"b_max\"" << "," <<

      "\"max_tuples_memory_size\"" << "," <<
      "\"tuples_processed\"" << "," <<
      "\"factor_calls\"" << "," <<
      "\"norm_equation_calls\"";

  return out.str();
}
