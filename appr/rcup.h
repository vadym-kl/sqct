#ifndef RCUP_H
#define RCUP_H

#include "tcount.h"
#include "hprhelpers.h"
#include "findhalves.h"

#include <array>

struct rcup
{
  typedef std::pair< hprr, min_unitaries > rcup_res;

  typedef std::pair< double, std::array<long,5> > res_tuple;
  typedef std::vector< res_tuple > res_tuples_arr;
  typedef std::pair<double,double> interval_t;


  rcup( long n, const hprr& phi, const hprr& delta );
  void merge_halves( const halves_t& re, const halves_t& im, int k );
  long recover_real( const halves_t::value_type& a, int k ) const;
  long recover_imag( const halves_t::value_type& a, int k ) const;

  rcup_res R;
  interval_t I;

  hprr cos2m[2];
  hprr sin2m[2];
  hprr reW[2];
  hprr imW[2];
  res_tuples_arr out;
};




#endif // RCUP_H
