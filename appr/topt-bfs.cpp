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

#include "topt-bfs.h"
#include "output.h"
#include "es/exactdecomposer.h"

#include "serializers.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>

#include <cassert>

using namespace std;


static void print_list( std::vector<int>& a, std::vector<int>& b, std::vector<int>& c)
{
    cout << "sde(|.|^2)   total   min-T   max-T" << endl;
    cout << "                     gates   gates" << endl;

    for( int i = 0; i < a.size(); ++i )
        if( a[i] != 0 )
            cout <<  setw (10) << i << ":" <<
                     setw (7) << a[i] << " " <<
                     setw (7) << b[i] << " " <<
                     setw (7) << c[i] << endl;
}

toptbfs::toptbfs() : og(max_cost)
{
    init();
}

//assumes that initial node is identity matrix
static void getCircuit( circuit& res, const optNode* node, const vector<int>& genIdToGateId )
{
    const optNode* current = node;
    while( current->parent != 0 )
    {
        res.push_front( genIdToGateId[current->gen_id] );
        current = current->parent;
    }
}

void toptbfs::init()
{
    static const gateLibrary& gl = gateLibrary::instance();
    const bool verbose  = true;
    typedef matrix2x2<int> m;

    ogc.m_generators = { m::X(), m::Y(), m::Z(), m::H(), m::P(), m::P().conjugateTranspose() };
    vector<int> genIdToGateId = { gl.X, gl.Y, gl.Z, gl.H, gl.P, gl.Pd };
    ogc.m_cost = { 1,2,1,10,40,40 };
    ogc.m_initial = { m::Id() };
    vector<int> initialGates = { gl.Id };
    ogc.m_initial_cost = {0};
    ogc.generate();

    m_clifford.resize( ogc.unique_elements().size() );
    vector<int> pauliCounts( ogc.unique_elements().size() );

    auto k = m_clifford.begin();
    auto k1 = pauliCounts.begin(); //number of Pauli gates used in each generator

    if( verbose ) cout << "Single qubit Clifford circuits:" << endl;

    for( auto i : ogc.unique_elements() )
    {
        getCircuit(*k,i.get(),genIdToGateId);
        auto counts = k->count();
        *k1 = counts[ gl.Z ] + counts[ gl.X ] + counts[ gl.Y ];

        if( verbose ) { k->toStreamSym(cout); cout << " " << *k1 << "," << endl; }

        og.m_generators.push_back( i->unitary * m::T() );
        og.m_cost.push_back( 1 );
        og.m_initial.push_back( i->unitary );
        ++k;++k1;
    }

    og.m_initial_cost.resize( og.m_initial.size(), 0 );
    og.generate();

    for( size_t k = 0; k < og.m_cost_stat.size(); ++k )
      if( og.m_cost_stat[k] > 0  )
        cout << k << ":" << og.m_cost_stat[k] << endl;
}

double dist( const pair<double,double>& x, double phi, int k )
{
  double precision = 1e-10;
  const double pi8 = M_PI / 8.0;
  double theta = pi8*k - phi*0.5;
  pair<double,double> w(cos(theta),sin(theta));
  double v1 = w.first * ( w.first - x.first ) + w.second * ( w.second - x.second );
  double v2 = w.first * ( w.first + x.first ) + w.second * ( w.second + x.second );
  if( v2 <= 0. && v2 >= -precision ) v2 = 0.;
  if( v1 <= 0. && v1 >= -precision ) v1 = 0.;
  return sqrt(min(v1,v2));
}

std::vector<std::pair< double, min_unitaries>> bfs_results::cup(const hprr &phi, int max_layer ) const
{
  std::vector<std::pair< double, min_unitaries> > R;
  R.resize(max_layer);
  double precision = 1e-10;

  R[0] =  cupl(phi,0,2.0);
  for( int layer = 1; layer < max_layer; ++layer )
    R[layer] = cupl(phi,layer,R[layer-1].first);

  return R;
}

std::pair< double, min_unitaries> bfs_results::cupl(const hprr &phi, int layer, double min_dist) const
{
  double precision = 1e-10;
  std::pair< double, min_unitaries> R;
  //finds all unitaries closest to R_z(phi) within a layer
  const float_index_t& ind = m_index[layer];
  const layer_t& lr = m_layers[layer];
  double md = min_dist;
  double ph = to_ld(phi);

  long id  = -1;
  for( size_t i = 0; i < ind.size(); ++i )
  {
//    cout << ind[i].x << "," << ind[i].k << endl;
//    cout << lr[ind[i].id] << endl;
    double d = dist(ind[i].x,ph,ind[i].k);
    if( d <= md + precision )
    {
      md = d;
      id = ind[i].id;
    }
  }

  R.first = md;

  if( id > 0 )
  {

    matrix2x2<int> m = global_phase_canonical(lr[id]);

    //find unitaries with the same 'x' entry in the layer
    ring_int<int> x = m.d[0][0];

    R.second.x = zwt(m.d[0][0]);
    R.second.m = m.de;
    R.second.k = m.det_power();

    int max_e = (1 << max_cost);
    ring_int<int> mx(max_e,max_e,max_e,max_e);
    int min_e = - max_e;
    ring_int<int> mn(min_e,min_e,min_e,min_e);

    for( int i = 0; i < 8; ++i )
    {
      matrix2x2<int> low(x,mn,mn,mn,R.second.m);
      matrix2x2<int> high(x,mx,mx,mx,R.second.m);
      auto lb = upper_bound(lr.begin(),lr.end(),low);
      auto ub = upper_bound(lr.begin(),lr.end(),high);

      for( ; lb != ub; ++lb )
      {
        matrix2x2<int> c = global_phase_canonical(*lb);
        if( c.d[0][0] == m.d[0][0] && c.det_power() == R.second.k )
          R.second.y.push_back(zwt(c.d[1][0].i_canonical()));
      }
      x.mul_eq_w();
    }

  }

  if( min_dist <= R.first + precision && min_dist >= R.first - precision )
  {
    R.second.x = zwt(0,0,0,0);
    R.second.y.clear();
  }
  else
  {
    R.second .min_t_count = layer;
    R.second.to_canonical_form();
  }
  return R;
}


string bfs_results::filename( int layer )
{
  stringstream ss;
  ss << "bfs-layer-" << layer;
  return ss.str();
}

const bfs_results &bfs_results::instance()
{
  static bfs_results res;
  ifstream ifs(filename(0)+".uni.bin");
  if( ifs.is_open() )
  {
    ifs.close();
    res.load();
  }
  else
  {
    res.get();
  }
  return res;
}


void bfs_results::get()
{
  m_layers.resize(max_cost);
  m_index.resize(max_cost);

  toptbfs tb;

  const auto& ue = tb.og.unique_elements();
  typedef optSequenceGenerator::nodeptr nd_t;
  for( const nd_t& a : ue )
  {
    if( a->cost < max_cost )
      m_layers[a->cost].push_back(global_phase_canonical(a->unitary));
  }

  for( int i = 0; i < max_cost; ++i )
  {
    std::sort( m_layers[i].begin(), m_layers[i].end() );
    write_stl( m_layers[i],filename(i) + ".uni.bin" );
  }

  m_index.resize(max_cost);
  for( int i = 0; i < max_cost; ++i )
  {
    const auto& lr = m_layers[i];
    auto& il = m_index[i];
    for( size_t k = 0; k < lr.size(); ++k )
    {
      matrix2x2hpr mf(lr[k]);
      il.push_back(float_index_entry{make_pair(to_ld(mf.d[0][0].real()),to_ld(mf.d[0][0].imag())),lr[k].det_power(),k});
    }
    write_stl( m_index[i],filename(i) + ".ind.bin" );
  }
}

void bfs_results::load()
{
  m_layers.resize(max_cost);
  m_index.resize(max_cost);
  for( int i = 0; i < max_cost; ++i )
  {
    read_stl( m_layers[i],filename(i) + ".uni.bin" );
    read_stl( m_index[i],filename(i) + ".ind.bin" );
  }
}
