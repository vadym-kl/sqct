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


#ifndef OBSERVATIONS_H
#define OBSERVATIONS_H

#include "hprhelpers.h"

#include <string>
#include <iostream>

struct observations
{
  observations();

  long n; //*
  hprr phi; //*
  hprr delta; //*

  double find_halves_real_time; //*
  double find_halves_cpu_time; //*
  double merge_halves_time; //*
  double tcount_time; //*

  size_t halves_size; //*
  long b_max;

  size_t max_tuples_memory_size; //*
  size_t tuples_processed; //*
  size_t factor_calls;
  size_t norm_equation_calls;

  static std::string title();

};

std::ostream& operator<< ( std::ostream& out, const observations& o );

#endif // OBSERVATIONS_H
