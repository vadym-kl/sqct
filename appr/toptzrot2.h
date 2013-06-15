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

#ifndef TOPTZROT2_H
#define TOPTZROT2_H

#include "approxlist.h"
#include "rint.h"
#include "matrix2x2.h"
#include "normsolver.h"
#include "symbolic_angle.h"
#include <tuple>
#include <string>


typedef matrix2x2<mpz_class> matr;

/// \brief Represents a unitary matrix
///  1 / sqrt{2}^n {{x,-y* w^k},{y, x* w^k}}, for w = Exp[i pi/4]
struct unitary
{
  unitary();
  explicit unitary( const matr& m );
  operator matr();
  ring_int<mpz_class> x;
  ring_int<mpz_class> y;
  int n;
  int k;
};

std::ostream& operator<<(std::ostream& out, const unitary& re );

/////////////////////////////////////////////////////////////////////

struct approx_params
{
  hprr phi;
  hprr epsilon;
  int sde2;
};

struct approx_result_block
{
  approx_result_block();
  std::vector<long double>  precisions;
  std::vector< matr >       matrices;
};

struct approx_result : public approx_result_block
{
  approx_result();
  void set_hint( int sde2, long double precision, const matr& m );
  void set_certified_epsilon( int sde2, long double eps );
  void set_result( int sde2, long double precision, const matr& m );

  approx_result_block hints; // ocasionally found results or results produced by other methods
  std::vector<long double>  certified_epsilon; //if certified_epsilon[k] > 0 there is no approximation better than epsilon below certified_epsilon[k] for sde2 = k
};

/// \brief Finds an approximation of R_z(<in.phi>) rotation by unitary over the ring Z[i,1/sqrt{2}].
///
/// Searches for unitaries with sde in the range [<sde2min>,<sde2max>].
/// Considers unitaries of the form 1/2^n {{x,-y* w^k},{y, x* w^k}}, for w = Exp[i pi/4] and x,y from Z[w].
/// Uses global phase invariant distance to determine the quality of the approximation.
///
/// This results in two approximation subproblems corresponding to
/// theta[0] = - <in.phi> / 2 and
/// theta[1] = - <in.phi> / 2 + pi / 8
///
/// Each subproblem is stated as following:
/// Find number x over the ring Z[i,1/sqrt{2}] such that:
///
struct toptzrot2
{
  toptzrot2(const approx_params& in, approx_result& out);
  void search();
  void combine_lists();
  long double search_initial_threshold();
  bool try_solution( long b, long d );
  bool done();
  void add_candidates_in_range( int i,long double eps_min, long double eps_max );

  const approx_params& in;
  approx_result& out;

  hprr        theta[2];
  approx_list cosTheta[2];///< cos(theta[i])
  approx_list sinTheta[2];///< sin(theta[i])
  hprr        delta;
  int         n;
  int         sde2min;
  int         sde2max;

  hprr        pow2nCosTh[2];
  hprr        pow2nSinTh[2];

  std::vector< std::tuple<long double, long,long > > candidates;

  std::vector<long double>  precisions;
  std::vector< matr >       matrices;
  matrix2x2hpr              target;
  const normSolver&                solver;
};

/////////////////////////////////////////////////////////////////////

struct optzrot_driver_params
{
  int sde2initial;
  int sde2max;
  hprr phi;
};

struct optzrot_driver
{
  optzrot_driver( const optzrot_driver_params& in, approx_result& out );

  const optzrot_driver_params& in;
  approx_result& out;
  void search();
};

struct result_entry
{
  result_entry();
  long double precision; //precision achieved by unitary
  long double cert_precision; //if there is no unitary, shows which precision was assumed
  bool cert_gap;
  long double required_precision;
  unitary u; //U[x,y,k]
  int t_gates; // number of gates we off from minimum for given x in U[x,y,k]
  int delta_t;
private:
  result_entry( const result_entry& );
};

std::ostream& operator<<(std::ostream& out, const result_entry& re );

struct topt_result : public std::vector< result_entry >
{
  topt_result();
  symbolic_angle angle;
  void write( const std::string& filename );
};

/// \brief tweaks results of toptzrot2 to achieve sde2 - 2 or sde2 - 1 T gates
struct toptimizer
{
  toptimizer( const approx_result& in, topt_result& out );
  unitary optimize_t_gates( const matr& m );
  const approx_result& in;
  topt_result&        out;
};

struct topt_app_params
{
  int max_sde;
  std::string in_filename;
  std::string out_filename;
};

struct topt_app
{
  topt_app( const topt_app_params& in ) :
    params(in) {}
  void run();
  const topt_app_params& params;
};

#endif // TOPTZROT2_H
