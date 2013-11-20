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


#ifndef RCUP_H
#define RCUP_H

#include "observations.h"

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
  interval_t IL;

  hprr cos2m[2];
  hprr sin2m[2];
  hprr reW[2];
  hprr imW[2];
  res_tuples_arr out;

  observations obs;
};




#endif // RCUP_H
