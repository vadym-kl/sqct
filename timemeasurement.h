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
  clock_gettime(CLOCK_REALTIME, &t);
}

#endif // TIMEMEASUREMENT_H
