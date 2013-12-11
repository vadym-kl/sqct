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

#ifndef TOPTBFS_H
#define TOPTBFS_H

#include "es/optsequencegenerator.h"
#include "gatelibrary.h"
#include "tcount.h"
#include "hprhelpers.h"

/// \brief Performes check of the conjecture of the T optimality of exact
/// decomposition of the algorithm
class toptbfs
{
public:
    static const int max_cost = 18;
    /// \brief Initialized by call of init()
    toptbfs() ;
    /// \brief Start exhaustive search of T-optimal circuits
    /// then applies exact synthesis algorithm to each of them
    /// and checks if the result of algorithm contains optimal
    /// number of T gates. Outputs to console number of non-optimal
    /// decompositions
    void init();
    /// \brief Optimal sequence generator used to get short circuits for Clifford group unitaries
    optSequenceGenerator ogc;
    /// \brief Optimal sequence generator used to get T-optimal circuits
    optSequenceGeneratorCostLim og;
    /// \brief Circuits for Clifford group unitaries
    std::vector< circuit > m_clifford;
};

struct float_index_entry
{
  std::pair<double,double> x;
  int k;
  long id;
};

struct bfs_results
{
  static const int max_cost = toptbfs::max_cost + 1;

  bfs_results() {}
  void get();
  void load();
  static std::string filename(int layer);

  static const bfs_results& instance();

  typedef std::vector< matrix2x2<int> > layer_t;
  typedef std::vector< float_index_entry > float_index_t;

  std::vector< layer_t > m_layers;
  std::vector< float_index_t > m_index;

  std::vector<std::pair< double, min_unitaries>> cup(const hprr& phi, int max_layer = max_cost ) const;
  std::pair< double, min_unitaries> cupl(const hprr& phi, int layer, double min_dist ) const;
};

#endif // TOPTBFS_H
