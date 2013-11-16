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

#include "toptzrot2.h"
#include <algorithm>
#include <fstream>
#include <cassert>
#include "output.h"
#include <tuple>

#include <ratio>
#include <chrono>

#include <omp.h>


using namespace std;
using namespace std::chrono;
typedef matrix2x2<mpz_class> matr;

//////////////////////////////////////////////////////////////

unitary::unitary(const matr &m) :
  x(m.d[0][0]),y(m.d[1][0]),n(m.de), k(0)
{
  matr sp(m);
  sp.d[1][1] = m.d[0][0].conjugate();
  sp.d[0][1] = - m.d[1][0].conjugate();

  for( int i = 0; i < 8; ++i )
  {
    if( sp * m.T(i) == m )
    {
      k = i;
      break;
    }
  }
}

unitary::operator matr()
{
  matr tmp(x, - y.conjugate() * y.omega(k) ,
           y,   x.conjugate() * y.omega(k) ,n);
  tmp.reduce();
  return tmp;
}

std::ostream& operator<<(std::ostream& out, const unitary& u )
{
  out << u.n << ",";
  for( int i = 0; i < 4; ++i )
    out << u.x[i] << ",";
  for( int i = 0; i < 4; ++i )
    out << u.y[i] << ",";
  out << u.k ;
  return out;
}

//////////////////////////////////////////////////////////////

static const int maximum_sde2 = 250;

struct time_logger
{
  time_logger( const string& event_name ) : name(event_name)
  {
    t1 = high_resolution_clock::now();
  }

  ~time_logger()
  {
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    ofstream ofs("time.log",ios_base::app );
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
    //ofs << "[" ; ofs.width(40);
    ofs << name << ",";
    //ofs.width(20);
    ofs << time_span.count() << endl;
  }

  high_resolution_clock::time_point t1;
  string name;
};

////////////////////////////////////

static long integer_part( const hprr& num, long b )
{
  hprr tmp = ( num - hprHelpers::sqrt2() * hprr(b) );
  if( ! mpfr_fits_slong_p( tmp._x, MPFR_RNDN ) )
    throw std::logic_error("overflow");
  return mpfr_get_si( tmp._x, MPFR_RNDN );
}

static ring_int_real<mpz_class> sqrt2_conjugate( const ring_int_real<mpz_class>& val )
{
  return ring_int_real<mpz_class>(val[0],-val[1]);
}

////////////////////////////////////

template< class T>
static inline void sort( vector< T >& vec )
{
  sort(begin(vec),end(vec));
}

static pair<long double,long double> vec_values_range( const vector< pair< long double, long> >& vec )
{
  return make_pair( vec.front().first, (vec.end() - 1)->first );
}

static inline auto lower_bound( const vector< pair< long double, long> >& vec, long double val ) -> decltype(vec.begin())
{
  return lower_bound(begin(vec),end(vec), make_pair(val,0l));
}

static inline auto upper_bound( const vector< pair< long double, long> >& vec, long double val ) -> decltype(vec.begin())
{
  return upper_bound(begin(vec),end(vec), make_pair(val,0l));
}

static long double vector_threshold( const vector< pair< long double, long> >& vec )
{
  pair<long double, long> zp(0,0);
  auto r1 = lower_bound( begin(vec), end(vec), zp );
  auto r2 = upper_bound( begin(vec), end(vec), zp );
  if( r1 == vec.end() || r2 == vec.end() )
    throw std::logic_error("initial threshold estimate failed");
  return max( abs(r1->first), abs(r2->first) );
}

///////////////////////////

toptzrot2::toptzrot2(const approx_params& _in, approx_result& _out) :
  in(_in), out(_out), precisions(maximum_sde2,2.0), matrices(maximum_sde2),
  solver( normSolver::instance() )
{
  theta[0] = - in.phi / hprHelpers::two();
  theta[1] = - in.phi / hprHelpers::two() + hprHelpers::pi() / hprr(8.0);
  n = ceil( (in.sde2 + 2.) / 4. );
  delta = in.epsilon * pow2(n) * sqrt(hprr(2.0));

  for( int i = 0; i < 2 ;++i )
  {
    pow2nSinTh[i] = pow2(n) * sin ( theta[i] );
    pow2nCosTh[i] = pow2(n) * cos ( theta[i] );
  }

  sde2min = 4 * n - 5;
  sde2max = 4 * n - 2;

  std::cout << "sde2 in [" << sde2min << "," << sde2max << "]" << endl;

  RzTh(theta[0],target);

  search();
}

void toptzrot2::search()
{

//before parallelization:
//  {
//    time_logger l("approx_list_builder," + to_string(n));
//    for( int i = 0; i < 2; ++i )
//    {
//      approx_list_builder lb0({n, cos(theta[i]),delta},cosTheta[i]);
//      approx_list_builder lb1({n, sin(theta[i]),delta},sinTheta[i]);
//    }
//  }


  {
    time_logger l("approx_list_builder," + to_string(n));
    int threads = omp_get_max_threads();
    cout << "Using " << threads << " threads" << endl;

    vector< approx_list > cl(threads * 2);
    vector< approx_list > sl(threads * 2);

    #pragma omp parallel
    //for( int thid = 0; thid < threads; ++thid )
    {
      int thid = omp_get_thread_num();
      for( int i = 0; i < 2; ++i )
      {
        approx_list_builder lb0({n, cos(theta[i]),delta,threads,thid},cl[ 2 * thid + i]);
        approx_list_builder lb1({n, sin(theta[i]),delta,threads,thid},sl[ 2 * thid + i]);
      }
    }

    size_t ssz[2] = {0,0};
    size_t csz[2] = {0,0};

    for( int th = 0; th < threads ; ++th )
    {
      for( int i = 0; i < 2; ++i )
      {
        csz[i] += cl[ 2 * th + i].size();
        ssz[i] += sl[ 2 * th + i].size();
      }
    }

    approx_list::iterator ci[2];
    approx_list::iterator si[2];

    for( int i = 0; i < 2; ++i )
    {
      cosTheta[i].resize(csz[i]);
      ci[i] = cosTheta[i].begin();

      sinTheta[i].resize(ssz[i]);
      si[i] = sinTheta[i].begin();
    }

    for( int th = 0; th < threads ; ++th )
    {
      for( int i = 0; i < 2; ++i )
      {
        ci[i] = copy(cl[ 2 * th + i].begin(),cl[ 2 * th + i].end(),ci[i]);
        si[i] = copy(sl[ 2 * th + i].begin(),sl[ 2 * th + i].end(),si[i]);
      }
    }

  }

  {
    time_logger l("sort," + to_string(n));
    for( int i = 0; i < 2; ++i )
    {
      sort(cosTheta[i]);
      sort(sinTheta[i]);
    }
  }

  {
    time_logger l("combine_lists," + to_string(n));
    combine_lists();
  }
}

void toptzrot2::combine_lists()
{
  bool not_done = true;
  long double search_threshold = to_ld( pow2(n) * in.epsilon * in.epsilon );

  long double eps1 = 0;
  long double eps2 = min( search_initial_threshold(), search_threshold );

  while( not_done && eps1 <= search_threshold )
  {
    candidates.clear();
    for( int i = 0; i < 2; ++i )
      add_candidates_in_range( i, eps1, eps2 );

    sort( candidates );
    for( const auto& t : candidates )
    {
      if( try_solution( get<1>(t), get<2>(t) ) )
        not_done = ! done();
    }

    eps1 = eps2;
    eps2 *= 2.0;
  }

  //if( !not_done )
  {
    cout << "epsilon = " << in.epsilon << endl;
    for( int k = sde2min; k <= sde2max; ++k )
    {
      cout << k << "," << precisions[k] ;
      auto cpr = precisions[k];
      if( cpr > in.epsilon )
      {
        if( cpr < 1. )
        {
          out.set_hint(k,precisions[k],matrices[k]);
          cout << "( > epsilon )" ;
        }
        out.set_certified_epsilon(k,to_ld(in.epsilon));
      }
      else
        out.set_result(k,precisions[k],matrices[k]);

      cout << endl;
    }

    for( size_t k = sde2max + 1; k < precisions.size(); ++k )
    {
      if( precisions[k] < 1. )
      {
        out.set_hint(k,precisions[k],matrices[k]);
        cout << "hint:" << k << "," << precisions[k] << endl;
      }
    }
    //return;
  }

  if( not_done )
    cout << "search threshold reached" << endl;
}

long double toptzrot2::search_initial_threshold()
{
  long double initial_threshold[2];
  for( int i = 0; i < 2; ++i )
  {
    initial_threshold[i] = max( vector_threshold(sinTheta[i]),
                                vector_threshold(cosTheta[i]));
  }
  return max( initial_threshold[0], initial_threshold[1] ) * 16.;
}

bool toptzrot2::try_solution( long b, long d )
{
  int theta_id = b & 1;
  //note: we use last bit of b and d to store which angle they approximate : Theta1 or Theta2
  hprr pow2nCT = pow2nCosTh[ theta_id ];
  hprr pow2nST = pow2nSinTh[ theta_id ];
  b >>= 1;
  d >>= 1;
  long a = integer_part( pow2nCT, b );
  long c = integer_part( pow2nST, d );

  mpz_class max = 1;
  max <<= 2 * n;

  ring_int<mpz_class> rr(a,d+b,c,d-b);
  auto rrabs2 = rr.abs2();

  if( rrabs2[0] > max )
    return false;

  if( abs(rrabs2.toComplex(0)) > pow2(2*n) )
    return false;

  if( abs(sqrt2_conjugate(rrabs2).toComplex(0)) > pow2(2*n) )
    return false;

  matr m;
  int sde2 = 4 * n - rrabs2.gde();

  if( sde2 < sde2min )
    return false;

  if(precisions[sde2] > 1. )
  {
    if( solver.solve(rr,n,m) )
    {
      if( theta_id == 0 )
        matrices[sde2] = m;
      else
        matrices[sde2] = m * m.T(1);


      precisions[sde2] = trace_dist(target,matrices[sde2]);
      return true;
    }
  }
  return false;
}

bool toptzrot2::done()
{
  for( int k = sde2min; k <= sde2max; ++k )
    if( precisions[k] > 1. )
      return false;
  return true;
}

void toptzrot2::add_candidates_in_range(int i, long double eps1, long double eps2)
{
  auto im_range = vec_values_range( sinTheta[i] );
  //      auto re_range = vec_values_range( cosTheta[i] );
  //      cout << im_range.first << ":" << im_range.second << endl;
  //      cout << re_range.first << ":" << re_range.second << endl;
  long double re_min = eps1 - im_range.second;
  long double re_max = eps2 - im_range.first;
  auto re_beg = lower_bound(cosTheta[i], re_min );
  auto re_end = upper_bound(cosTheta[i], re_max );
  for( auto k = re_beg; k < re_end; ++k )
  {
    long double im_min = eps1 - k->first;
    long double im_max = eps2 - k->first;
    auto im_beg = lower_bound( sinTheta[i], im_min );
    auto im_end = upper_bound( sinTheta[i], im_max );
    for( auto j = im_beg; j < im_end; ++j )
      candidates.push_back( make_tuple( k->first + j->first,
                                        (k->second << 1) + i,
                                        (j->second << 1) + i  ) );
  }
}



approx_result_block::approx_result_block() :
  precisions(maximum_sde2,2.),
  matrices(maximum_sde2)
{
}


approx_result::approx_result() :
  certified_epsilon(maximum_sde2,2.)
{
}

void approx_result::set_hint(int sde2, long double precision, const matr &m)
{
  long double cpr = hints.precisions[sde2];
  if( cpr < precision )
  {
    hints.precisions[sde2] = precision;
    hints.matrices[sde2] = m;
  }
}

void approx_result::set_certified_epsilon(int sde2, long double eps)
{
  if( certified_epsilon[sde2] > eps )
    certified_epsilon[sde2] = eps;
}

void approx_result::set_result(int sde2, long double precision, const matr &m)
{
  precisions[sde2] = precision;
  matrices[sde2] = m;
}


optzrot_driver::optzrot_driver(const optzrot_driver_params &_in, approx_result &_out) :
  in(_in), out(_out)
{
  search();
}

static void print_vec( const std::vector<long double>& v)
{
  for( auto i : v )
    cout << i << "," ;
  cout << endl;
}

static long double update( long double epsilon, int n )
{
  long double delta = to_ld( hprr(epsilon) * pow2(n) * hprHelpers::sqrt2() );
  if( delta > 0.5 )
  {
    epsilon = to_ld( hprr(0.4999999) / ( pow2(n) * hprHelpers::sqrt2() ) );
    cout << "epsilon was updated to:" << epsilon << endl;
  }
  return epsilon;
}

void optzrot_driver::search()
{
  out.precisions[0] = trace_dist(matrix2x2hpr::Id(), RzTh(in.phi) );
  out.matrices[0] = matr::Id();

  int n = ceil( ( in.sde2initial + 2.0 ) / 4. );
  int n_max = ceil( ( in.sde2max + 2.0 ) / 4. );

  for( ; n <= n_max; ++n )
  {
    long double epsilon = *min_element(out.precisions.begin(),
                          out.precisions.begin() + 4*n - 1 );
    cout << "epsilon:" << epsilon << endl;
    epsilon = update( epsilon, n ); //in the case if algorithm is not applicable for given precision

    approx_params inp{in.phi,epsilon, 4 * n - 2 };
    try
    {

      toptzrot2 tr(inp,out);
    }
    catch(const std::exception& e )
    {
      cerr << in.phi << "," << n << "," << epsilon << "," << ":" << e.what() << endl;
    }
    //here we need to decide which epsilon to feed next

    cout << "epsilon:" << epsilon << endl;
  }

}


unitary::unitary() :
  x(1,0,0,0), y(0,0,0,0), n(0), k(0)
{

}



toptimizer::toptimizer(const approx_result &_in, topt_result& _out ) :
  in(_in), out(_out)
{
  for (size_t i = 0; i < out.size(); ++i)
  {
    auto& o = out[i];

    auto eps = min_element(out.begin(), out.begin() + i,
                [] ( const result_entry& a, const result_entry& b )
                {
                  return a.precision < b.precision;
                }
               )->precision;

    if( in.precisions[i] < eps )
    {
      o.precision = in.precisions[i];
      o.u = optimize_t_gates(in.matrices[i]);
      matr m(o.u);
      o.t_gates = m.t();
      o.delta_t = m.t() - m.h() + 1;
    }
    else if( in.precisions[i] < 1.0 )
    {
      o.cert_precision = eps;
      o.cert_gap = false;
    } else {
      if( in.certified_epsilon[i] < 1.0 )
      {
        o.cert_precision = in.certified_epsilon[i];
        o.cert_gap = in.certified_epsilon[i] < eps ;
        if( o.cert_gap )
          o.required_precision = eps;
      }
    }
  }
}

unitary toptimizer::optimize_t_gates( const matr& m )
{
  int min_i = 0;
  bool is_conjugate = false;
  unitary u(m);

  int t_min = m.t();

  for( int i = 0; i < 8; ++i )
  {
    u.y.mul_eq_w();
    int tc = matr(u).t();
    if( t_min > tc )
    {

      t_min = tc;
      min_i = i + 1;
    }
  }

  u.y.conjugate_eq();

  for( int i = 0; i < 8; ++i )
  {
    u.y.mul_eq_w();
    int tc = matr(u).t();
    if( t_min > tc )
    {
      t_min = tc;
      min_i = i + 1;
      is_conjugate = true;
    }
  }

  u = (unitary)m;
  if( is_conjugate )
    u.y.conjugate_eq();
  u.y.mul_eq_w(min_i);

  return u;
}


topt_result::topt_result() :
  vector<result_entry>(maximum_sde2)
{
}

void topt_result::write(const string &filename)
{
  size_t max_items = 0;
  for( size_t i = 0; i < size(); ++i )
    if( at(i).precision < 1.0 )
      max_items = i + 1;

  ofstream fout(filename + ".data", ios_base::app );
  ofstream fc(filename , ios_base::app );
  ofstream ftc(filename + ".tc.csv", ios_base::app );

  for( size_t i = 0; i < max_items; ++i )
  {
    if( at(i).precision < 1.0 || at(i).cert_gap )
    {
      fout << angle << "," << i << "," << at(i) << endl;
      if( at(i).precision < 1.0 )
      {
        fc << angle.numerator << "," << angle.denominator << "," <<
              at(i).precision << "," << exactDecomposer::decompose(at(i).u) << endl;
      }

      if( at(i).precision < 1.0 && at(i).cert_gap )
      {
          ftc <<  at(i).precision << "," << at(i).t_gates << endl;
      }
    }
  }
}





result_entry::result_entry() :
  precision(2.0),
  cert_precision(2.0),
  cert_gap(false),
  required_precision(2.0),
  t_gates(0),
  delta_t(0)
{
}

std::ostream& operator<<(std::ostream& out, const result_entry& re )
{
  out << re.precision << "," << re.t_gates << "," << re.delta_t << ","
      << re.cert_precision << "," << re.cert_gap << "," << re.required_precision << ","
      << re.u;
  return out;
}





void topt_app::run()
{
  ifstream ifs;
  ifs.exceptions (  std::ifstream::badbit );
  ifs.open(params.in_filename);
  vector<symbolic_angle> angles;
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
    ss >> a.numerator >> a.denominator >> a.pi;
    angles.push_back(a);
  }

  cout << "scheduled to find approximations for " << angles.size() << " unitaries";

  for( const symbolic_angle& a : angles )
  {
    optzrot_driver_params ap;

    ap.sde2initial = 15;
    ap.sde2max = params.max_sde;
    ap.phi = a;

    approx_result ares;
    optzrot_driver drv(ap,ares);

    topt_result out;
    out.angle = a;
    toptimizer tpz(ares,out);
    out.write(params.out_filename);
  }
}
