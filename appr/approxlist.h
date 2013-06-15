//     Copyright (c) 2013 Vadym Kliuchnikov sqct(dot)software(at)gmail(dot)com
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
//     Based on "Practical approximation of single-qubit unitaries by single-qubit quantum Clifford and T circuits"
//     by Vadym Kliuchnikov, Dmitri Maslov, Michele Mosca
//     http://arxiv.org/abs/1212.6964v1, new version of the paper describing this modification: [to be published]

#ifndef APPROXLIST_H
#define APPROXLIST_H

#include "hprhelpers.h"
#include <vector>

/// \brief Specifies input for approx_list_builder struct.
/// See approx_list_builder struct for the meaning of each member.
struct approxlist_params
{
  int n;
  hprr s;
  hprr delta;
  int step;
  int init;

  approxlist_params( int _n , hprr _s, hprr _delta, int _step = 1, int _init = 0 ) :
    n( _n ), s(_s), delta(_delta), step(_step), init(_init)
  {}

};

typedef std::vector< std::pair< long double, long > > approx_list;

/// \brief Finds an approximations of a real number <in.s>
/// using numbers of the from ( a + \sqrt{2} b ) / 2^<in.n>;
/// b is chosen from the set  ( b_min + <in.init> + N * <in.step> ), where N is a non-negative integer;
/// b_min calculated based on <in.n> and <in.s>;
/// a is an integer;
///
/// <in.init> and <in.step> members we introduced to run approximation using multiple threads
///
/// For any approximation found it is true that:
/// | a + \sqrt{2} b - 2^<in.n> * s | <= <in.delta>
/// and a^2 + b^2 <= 4^<in.n>
///
/// Implemented algorithm supports only <in.delta> < 0.5. In this case a uniquely defined by b
///
/// The result of the execution is a vector of pairs:
/// ( ( 2^<in.n> * <in.s>  - ( a + \sqrt{2} b ) ) * <in.s> , b )
///
/// in.<name> refers to the member named <name> of approxlist_params struct
struct approx_list_builder
{
  approx_list_builder( const approxlist_params& in,
                       approx_list& out );
  bool norm_condition( long a, long b );

  const approxlist_params& in;
  approx_list& out;
  __int128 pow4n;

};

#endif // APPROXLIST_H
