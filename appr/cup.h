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


#ifndef CUP_H
#define CUP_H

#include "hprhelpers.h"
#include "tcount.h"

#include <string>



struct cup_params
{
  static const int default_max_lookup = 8;
  const hprr phi;
  int max_layer;
  int max_lookup;
  std::string angle_str;
  std::string results_file;
};

struct cup
{
public:
  typedef std::pair<double,min_unitaries> result_t;
  cup(const hprr &phi, int max_layer, int max_lookup);
  cup( const cup_params& params);
  std::vector< result_t > R;
};

#endif // CUP_H
