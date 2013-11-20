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

#ifndef TIMEMEASUREMENT_H
#define TIMEMEASUREMENT_H

#include <time.h>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>


typedef boost::multiprecision::int512_t i512;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<200> > f200;
typedef timespec timepoint;

static i512 time_diff_i( const timepoint& t1, const timepoint& t2)
{
  i512 nano_inv = 1000000000;
  i512 T1 = i512(t1.tv_sec) * nano_inv + i512(t1.tv_nsec);
  i512 T2 = i512(t2.tv_sec) * nano_inv + i512(t2.tv_nsec);
  return T2 - T1;
}

static double time_diff_d( const timepoint& t1, const timepoint& t2)
{
  i512 dt = time_diff_i(t1,t2);
  i512 k = 1000000; //get miliseconds
  f200 r = f200(dt) / f200(k);
  return r.convert_to<double>();
}

static void get_timepoint( timepoint& t )
{
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t);
}

static void get_timepoint_real( timepoint& t )
{
  clock_gettime(CLOCK_REALTIME, &t);
}

#endif // TIMEMEASUREMENT_H
