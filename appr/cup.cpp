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

#include "cup.h"
#include "topt-bfs.h"

#include "rcup.h"
#include "matrix2x2.h"

#include "es/exactdecomposer.h"

#include <iostream>
#include <fstream>

using namespace std;

cup::cup(const hprr &phi, int max_layer, int max_lookup)  :
  R(max_layer)
{
  const bfs_results& br = bfs_results::instance();

  auto r = br.cup(phi,max_lookup);
  const int mc = bfs_results::max_cost;
  int bnd = std::min(max_layer,std::min(max_lookup,mc));
  for( int i = 0; i < bnd ; ++i )
  {
    R[i] = r[i];
    if( ! R[i].second.y.empty() )
      R[i].first = trace_dist(Rz(phi),matrix2x2<mpz_class>(R[i].second));
  }

  bool verbose = false;

  if( verbose )
    cout << "{" << observations::title() << "}," << endl;

  for( int i = bnd; i < max_layer; ++i )
  {
    rcup rc(i,phi,R[i-1].first);
    if( verbose )
      cout << "{" << rc.obs << "}," << endl;
    R[i] = std::make_pair(to_ld(rc.R.first),rc.R.second);

    if( ! R[i].second.y.empty() )
      R[i].first = trace_dist(Rz(phi),matrix2x2<mpz_class>(R[i].second));
  }
}

std::string circuit_title()
{
  stringstream ss;
  ss << "\"Matrix product\",\"T-count\",\"H-count\",\"Phase-count\",\"Pauli-X-count\",\"Pauli-Y-count\",\"Pauli-Z-count\"" ;
  return ss.str();
}

std::string ciruit_str( const matrix2x2<mpz_class>& m )
{
  const exactDecomposer& es = exactDecomposer::instance();
  circuit res = es.decompose(m);
  auto cost = gateLibrary::toCliffordT(res.count());
  stringstream ss;

  ss << "\"";
  res.toMathStream(ss);
  ss << "\"," << cost[gateLibrary::T] << ","
     << cost[gateLibrary::H] << ","
     << cost[gateLibrary::P] << ","
     << cost[gateLibrary::X] << ","
     << cost[gateLibrary::Y] << ","
     << cost[gateLibrary::Z] ;
  return ss.str();
}

std::string ciruit_str_empty()
{
  stringstream ss;

  ss << "Id,0,"
     << "0" << ","
     << "0" << ","
     << "0" << ","
     << "0" << ","
     << "0" ;
  return ss.str();
}

void write_title( const cup_params& params, ostream& ofs )
{
  ofs << "{";
  ofs << "\"Angle string\","
      << "\"Angle value\","
      << "\"T count\","
      << "\"Precision\","
      << "{" << min_unitaries::short_title() << "},"
      << circuit_title() << ",";
  ofs << observations::title() << "}" << endl;
}

void write_line( int t_count, const cup_params& params, ostream& ofs, const typename cup::result_t& res, observations& obs )
{
  ofs << "{";
  ofs << params.angle_str << ",";
  ofs << fixed;
  ofs.precision(30);
  ofs << to_ld(params.phi) << "," <<
         t_count << "," <<
         res.first << "," <<
         "{" << res.second.short_str() << "},";
  if( res.second.y.empty() )
     ofs << ciruit_str_empty() << ",";
  else
     ofs << ciruit_str(res.second) << ",";
  ofs << obs << "}" << endl;
}

cup::cup(const cup_params &params)  :
  R(params.max_layer)
{
  const hprr &phi = params.phi;
  int max_layer = params.max_layer;
  int max_lookup = params.max_lookup;

  const bfs_results& br = bfs_results::instance();

  ofstream ofs(params.results_file,ios_base::app);

  {
    ofstream ofst(params.results_file + ".title");
    write_title(params,ofst);
  }

  auto r = br.cup(phi,max_lookup);
  const int mc = bfs_results::max_cost;
  int bnd = std::min(max_layer,std::min(max_lookup,mc));

  observations empty;
  for( int i = 0; i < bnd ; ++i )
  {
    R[i] = r[i];
    if( ! R[i].second.y.empty() )
    {
      R[i].first = trace_dist(Rz(phi),matrix2x2<mpz_class>(R[i].second));
      R[i].second.to_canonical_form();
    }
    write_line(i,params,ofs,R[i],empty);
  }

  for( int i = bnd; i < max_layer; ++i )
  {
    try
    {
      rcup rc(i,phi,R[i-1].first);
      R[i] = std::make_pair(to_ld(rc.R.first),rc.R.second);

      if( ! R[i].second.y.empty() )
      {
        R[i].first = trace_dist(Rz(phi),matrix2x2<mpz_class>(R[i].second));
        R[i].second.to_canonical_form();
      }

      write_line(i,params,ofs,R[i],rc.obs);
    }
    catch(...)
    {
      cout << "something went wrong:" << i << "," << phi << "," << params.angle_str <<endl;
    }
  }

}
