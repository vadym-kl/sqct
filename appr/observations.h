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

  double find_halves_time; //*
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
