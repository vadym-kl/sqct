#ifndef CUP_H
#define CUP_H

#include "hprhelpers.h"
#include "tcount.h"

struct cup
{
public:
  typedef std::pair<double,min_unitaries> result_t;
  cup( const hprr& phi, int max_layer, int max_lookup = 13 );
  std::vector< result_t > R;
};

#endif // CUP_H
