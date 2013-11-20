#include "cup.h"
#include "topt-bfs.h"

#include "rcup.h"

#include <iostream>

using namespace std;

cup::cup(const hprr &phi, int max_layer, int max_lookup ) :
  R(max_layer)
{
  const bfs_results& br = bfs_results::instance();

  auto r = br.cup(phi,max_lookup);
  const int mc = bfs_results::max_cost;
  int bnd = std::min(max_layer,std::min(max_lookup,mc));
  for( int i = 0; i < bnd ; ++i )
    R[i] = r[i];

  //TODO: add distance reevaluation -- to ensure accuracy

  bool verbose = true;

  if( verbose )
    cout << "{" << observations::title() << "}," << endl;

  for( int i = bnd; i < max_layer; ++i )
  {
    rcup rc(i,phi,R[i-1].first);
    if( verbose )
      cout << "{" << rc.obs << "}," << endl;
    R[i] = std::make_pair(to_ld(rc.R.first),rc.R.second);
  }
}
