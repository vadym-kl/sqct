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



#include "observations.h"

#include <sstream>

using namespace std;

std::ostream &operator <<(std::ostream &out, const observations &o)
{
  out << std::fixed;
  out.precision(30);

  out << o.n << "," <<
      to_ld(o.phi) << "," <<
      to_ld(o.delta) << "," <<

      o.find_halves_real_time << "," <<
      o.find_halves_cpu_time << "," <<
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

  find_halves_real_time(0),
  find_halves_cpu_time(0),
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

      "\"find_halves_real_time\"" << "," <<
      "\"find_halves_cpu_time\"" << "," <<
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
