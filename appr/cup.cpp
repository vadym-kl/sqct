#include "cup.h"
#include "topt-bfs.h"

#include "rcup.h"

cup::cup(const hprr &phi, int max_layer, int max_lookup ) :
  R(max_layer)
{
  bfs_results br;
  br.load();
  auto r = br.cup(phi);
  const int mc = bfs_results::max_cost;
  int bnd = std::min(max_layer,std::min(max_lookup,mc));
  for( int i = 0; i < bnd ; ++i )
    R[i] = r[i];

  //TODO: add distance reevaluation -- to ensure accuracy

  for( int i = bnd; i < max_layer; ++i )
  {
    rcup rc(i,phi,R[i-1].first);
    R[i] = std::make_pair(to_ld(rc.R.first),rc.R.second);
  }
}
