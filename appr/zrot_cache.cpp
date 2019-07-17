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

#include "zrot_cache.h"

#include <fstream>
#include <sstream>
#include <algorithm>
#include "output.h"

using namespace std;


static void normalize( step_function &sf1 )
{
  step_function tmp;
  vector<int> tc(sf1.size());
  for( size_t k = 0; k < sf1.size(); ++k )
    tc[k] = ((matr)sf1[k].second).t();

  tmp.push_back(sf1[0]);
  for( size_t k = 1; k < sf1.size(); ++k )
  {
    if( tc[k - 1] > tc[k] )
      tmp.push_back(sf1[k]);
  }

  swap(tmp,sf1);
}

zrot_cache::zrot_cache() : recalculate(false)
{}

zrot_cache::zrot_cache(const string &filename, bool rec ) : recalculate(rec)
{
  read_file(filename);
}

void zrot_cache::read_file(const std::string &filename)
{
  ifstream ifs;
  ifs.exceptions (  std::ifstream::badbit );
  ifs.open(filename);
  while( ifs )
  {
    string line;
    getline(ifs,line);

    if(line.size() == 0 )
      continue;
    if( line[0] == '#' )
      continue;

    //cout << line << endl;

    istringstream ss(line);
    symbolic_angle a{0,1,false};
    circuit c;
    double file_precision;
    char delim;
    ss >> a.numerator >> delim >> a.denominator >> delim >> file_precision >> delim;
    c.fromStream(ss,false);
    if(recalculate)
      (*this)[a].push_back({trace_dist(Rz(a),c),c});
    else
      (*this)[a].push_back({file_precision,c});
  }

  for ( auto& sf : *this )
  {
    sort(sf.second.begin(),sf.second.end());
    normalize(sf.second);
  }
}

circuit zrot_cache::lookup(const symbolic_angle &angle, double precision) const 
{
  auto it = find(angle);
  if ( it == end() )
    throw not_found_exception(angle,1.);

  const auto& sf = it->second;
  step_function::value_type v{precision,{}};
  auto lb = std::lower_bound( std::begin(sf), std::end(sf), v );
  if( lb == std::begin(sf) )
    throw not_found_exception(angle,precision);

  return (lb - 1)->second;
}

std::ostream& operator<< ( std::ostream& out , const zrot_cache& zc )
{
  out << "Precision, T count, circuit" << endl;
  for( const auto& v : zc )
  {
    out << v.first << ":" << endl;

    for( const auto& r : v.second )
    {
        int TC = gateLibrary::toCliffordT(r.second.count())[gateLibrary::T];
        //out << "    " <<  r.first << "," << TC << "," << r.second << endl;
        out.setf ( ios_base::fixed );
        out.precision( 100 );
        out << "{" <<  r.first << "," << TC << "},"<< endl;
    }
  }

  return out;
}




void DiffApplication::run()
{
  if( params.filenames.size() != 2 )
    throw logic_error("wrong number of inputs");
  zrot_cache zc1(params.filenames[0],true);
  zrot_cache zc2(params.filenames[1],true);

  cout << "File contents:" << endl;

  cout << "file:" << params.filenames[0] << endl;
  cout << zc1 << endl;
  cout << "**************************" << endl;

  cout << "file:" << params.filenames[1] << endl;
  cout << zc2 << endl;
  cout << "**************************" << endl;

  // find the angles that are common for both files

  vector<symbolic_angle> angles;

  for( const auto& rec : zc1 )
  {
    if(zc2.count(rec.first) > 0 )
      angles.push_back(rec.first);
  }

  cout << "Angles in common:" << angles.size() << endl;
  for( const auto& angle : angles )
    cout << "  " << angle << endl;

  for( const auto& angle : angles )
    compare(zc1[angle],zc2[angle],angle);

}

const long double low_k =  0.99999999;
const long double high_k = 1.00000001;



static bool belongs( const step_function &sf1, long double val)
{
  auto en = upper_bound(sf1.begin(), sf1.end(), pair<long double,circuit>{val * low_k,{}} );
  if( en == sf1.end() ) return false;
  return en->first < val * high_k;
}

static const circuit& circ( const step_function &sf1, long double val )
{
  auto en = upper_bound(sf1.begin(), sf1.end(), pair<long double,circuit>{val * low_k,{}} );
  return en->second;
}


static bool compare_circuits( const circuit& a, const circuit& b )
{
  matr ma(a);
  matr mb(b);
  ma.reduce();
  mb.reduce();

  if( ma.h() != mb.h() )
  {
    cout << "totally wrong: different sde|.|^2" << endl;
    throw std::logic_error("");
    return false;
  }

  bool found = false;
  int glob_ph = 0;
  for( int i = 0; i < 8; ++i )
  {
    if ( ma.d[0][0] == mb.d[0][0] * mb.d[0][0].omega(i) )
    {
      found = true;
      glob_ph = i;
    }
  }

  if( !found )
  {
    cout << "unexpected: different x" << endl;
    throw std::logic_error("");
    return false;
  }

  mb.mul_eq_w(glob_ph);

  if( mb.t() != ma.t() )
  {
    cout << endl << "different number of t gates:" << ma.t() << " vs " << mb.t() << endl;
    if(abs( mb.t() - ma.t() ) > 1 )
      throw logic_error("too big difference in the number of T gates");
    return false;
  }

  return true;
}

//template< class T>
//ostream& operator<< ( ostream& out, const vector<T>& vec )
//{
//  for( const auto& a : vec )
//    out << a << " ,";
//  return out;
//}


template< class T>
ostream& operator<< ( ostream& out, const vector< pair<T,circuit> >& vec )
{
  for( const auto& a : vec )
    out << a.first << " ,";
  return out;
}


void DiffApplication::compare(const step_function &sf1, const step_function &sf2, const symbolic_angle &angle)
{

  double dist_id = trace_dist(Rz(angle),matrix2x2hpr::Id()) * low_k;
  cout << "comparing angle:" << angle << " " ;
  pair<long double,circuit> vid {dist_id,{}};
  long double max1  = ( lower_bound( sf1.begin(), sf1.end(), vid ) - 1 )->first ;
  long double max2 = ( lower_bound( sf2.begin(), sf2.end(), vid ) - 1 )->first ;

  long double minl = max( sf1.begin()->first, sf2.begin()->first ) * low_k;
  long double maxl = min(max1,max2) * high_k ;

  auto bg = lower_bound(sf1.begin(), sf1.end(), pair<long double,circuit>{minl,{}} );
  auto en = upper_bound(sf1.begin(), sf1.end(), pair<long double,circuit>{maxl,{}} );

  std::vector<long double> mfs;
  std::vector<long double> diff_circs;
  int records = 0;
  for( auto i = bg; i < en; ++i  )
  {
    if( belongs(sf2,i->first) )
    {
      records++;
      if( !compare_circuits(i->second,circ(sf2,i->first)) )
        diff_circs.push_back(i->first);
    }
    else
    {
      mfs.push_back(i->first);
    }
  }

  bg = lower_bound(sf2.begin(), sf2.end(), pair<long double,circuit>{minl,{}} );
  en = upper_bound(sf2.begin(), sf2.end(), pair<long double,circuit>{maxl,{}} );

  std::vector<long double> mff;
  for( auto i = bg; i < en; ++i  )
  {
    if( !belongs(sf1,i->first) )
    {
      mff.push_back(i->first);
    }
  }

  cout << " total:" << records << " ";
  cout << "different in t gates: " <<  diff_circs.size() ;
  cout << " missing: from 1st: " << mff.size() << " ,from 2nd: " << mfs.size() << endl;

  if( mff.size() > 0 )
  {
    cout << "   missed from 1st: " << mff << endl;
  }

  if( mfs.size() > 0 )
  {
    cout << "   missed from 2nd: " << mfs << endl;
  }



  }


  void sqct_light_app::run()
  {
    ifstream ifs;
    ifs.exceptions (  std::ifstream::badbit );
    ifs.open(par.in_filename);
    vector<symbolic_angle> angles;
    vector<double> precisions;
    vector<string> names;
    while( ifs )
    {
      string line;
      getline(ifs,line);
      if(line.size() == 0 )
        continue;
      if( line[0] == '#' )
        continue;

      istringstream ss(line);
      symbolic_angle a{0,1,true};
      //char delim;
      double precision = 0.;
      string name;
      ss >> a.numerator >> a.denominator >> a.pi >> precision >> name;
      angles.push_back(a);
      precisions.push_back(precision);
      names.push_back(name);
    }

    cout << "Loading cahce ... ";
    zrot_cache zc(par.cache_filename);

    cout << "Processing " << angles.size() << " requests" << endl;

    ofstream ofs;
    ifs.exceptions (  std::ifstream::badbit );
    ofs.open(par.out_filename);


    for( size_t i = 0; i < angles.size(); ++i )
    {
      try
      {
        circuit c = zc.lookup(angles[i],precisions[i]);
      }
      catch(not_found_exception& e)
      {
        cout << "Not found:" << e.angle << ":" << e.precision << endl;
      }
    }

    //write footer
  }
