#ifndef FINDHALVES_H
#define FINDHALVES_H

#include "hprhelpers.h"

#include <vector>

hprr weight( const hprr& alpha, int m );

typedef std::vector< std::pair< double, long > > halves_t;

halves_t findhalves( const hprr& alpha, int m, const hprr& delta );

halves_t findhalves2( const hprr& alpha, int m, const hprr& delta );

halves_t findhalves3( const hprr& alpha, int m, const hprr& delta );

#endif // FINDHALVES_H
