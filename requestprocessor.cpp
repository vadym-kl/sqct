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


#include "requestprocessor.h"
#include "hprhelpers.h"
#include "appr/cup.h"

#include <iostream>
#include <sstream>
#include <map>

using namespace std;


typedef void request_processor( const vector<string>& request );

void process_pow2_request( const vector<string>& request )
{
  if( request.size() < 6 )
  {
    cout << "not enough parameters" << endl;
    return;
  }
  cout << "Scheduled approximation of R_z(Pi/2^k)" << endl;

  int Tcount_min = 0;
  int Tcount_max = 80;

  long k_min = 0;
  long k_max = 10;

  string filename = request[1];

  istringstream(request[2]) >> Tcount_min;
  istringstream(request[3]) >> Tcount_max;

  istringstream(request[4]) >> k_min;
  istringstream(request[5]) >> k_max;

  for( long k = k_min; k < k_max; ++k )
  {
    hprr phi = hprHelpers::pi() / pow2(k);
    stringstream astr;
    astr << "\"Pi/(2^" << k << ")\"" ;
    cup_params cp{phi,Tcount_max,cup_params::default_max_lookup,astr.str(),filename};
    cup c(cp);
  }
}

void process_uniform_request( const vector<string>& request )
{
  if( request.size() < 7 )
  {
    cout << "not enough parameters" << endl;
    return;
  }
  cout << "Scheduled approximation of R_z(2*Pi*k/n)" << endl;

  int Tcount_min = 0;
  int Tcount_max = 80;

  long n = 0;
  long k_min = 0;
  long k_max = 10;

  string filename = request[1];

  istringstream(request[2]) >> Tcount_min;
  istringstream(request[3]) >> Tcount_max;

  istringstream(request[4]) >> n;

  istringstream(request[5]) >> k_min;
  istringstream(request[6]) >> k_max;

  for( long k = k_min; k < k_max; ++k )
  {
    hprr phi = hprr(2)*hprHelpers::pi()*hprr(k) / hprr(n);
    stringstream astr;
    astr << "\"2*Pi*" << k << "/" << n << "\"";
    cup_params cp{phi,Tcount_max,cup_params::default_max_lookup,astr.str(),filename};
    cup c(cp);
  }
}

void process_angles_request( const vector<string>& request )
{
  cout << "Scheduled approximation R_z(alpha)" << endl;

  if( request.size() < 4 + 1)
  {
    cout << "not enough parameters" << endl;
    return;
  }

  int Tcount_min = 0;
  int Tcount_max = 80;

  string filename = request[1];

  istringstream(request[2]) >> Tcount_min;
  istringstream(request[3]) >> Tcount_max;

  for( size_t k = 4; k < request.size(); ++k )
  {
    hprr angle ;
    string angle_string;
    istringstream(request[k]) >> angle >> angle_string;
    cup_params cp{angle,Tcount_max,cup_params::default_max_lookup,angle_string,filename};
    cup c(cp);
  }
}
////////////////////////////////////////////////////////////////////

int process_request( const vector<string>& request )
{
  if(request.size() == 0)
  {
    cout << "Empty request" << endl;
    return 0;
  }
  map<string,request_processor*> reques_types_map;

  reques_types_map["ANGLES"] = process_angles_request;
  reques_types_map["UNIFORM"] = process_uniform_request;
  reques_types_map["POW2"] = process_pow2_request;

  if( reques_types_map.count(request[0]) != 0 )
    reques_types_map[request[0]](request);
  else
    cout << "unsupported request type:" << request[0] << endl;

}
