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

#ifndef ZROT_CACHE_H
#define ZROT_CACHE_H

#include "toptzrot2.h"
#include "gatelibrary.h"
#include <string>
#include <vector>
#include <stdexcept>

struct not_found_exception : public std::logic_error
{
  not_found_exception( const symbolic_angle& _angle, double _precision )
    : logic_error("approximation not found"),
      angle(_angle),
      precision(_precision) {};

  symbolic_angle angle;
  double precision;
};

typedef std::vector<std::pair<long double,circuit> > step_function;

struct zrot_cache : public std::map< symbolic_angle, step_function >
{
  zrot_cache();
  zrot_cache( const std::string& filename, bool recalculate = false );
  // loads cached data from file
  void read_file( const std::string& filename );
  // looks up approximation
  circuit lookup(const symbolic_angle& angle, double precision ) const ;
  bool recalculate;
};

std::ostream& operator<< ( std::ostream& out , const zrot_cache& zc );

struct DiffApplicationParams
{
  std::vector< std::string  > filenames;
};

struct DiffApplication
{
  DiffApplication( const DiffApplicationParams& in ) : params(in) {}
  void run();
  void compare(const step_function& sf1, const step_function& sf2, const symbolic_angle& angle );

  const DiffApplicationParams& params;
};

struct sqct_light_params
{
  std::string in_filename;
  std::string out_filename;
  std::string cache_filename;
};

struct sqct_light_app
{
  sqct_light_app( const sqct_light_params& in ) :
    par(in) {}
  void run();
  const sqct_light_params& par;
};

#endif // ZROT_CACHE_H
