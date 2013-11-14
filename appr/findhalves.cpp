
#include "findhalves.h"
#include <algorithm>

typedef hprr ft;

using namespace std;

hprr weight( const hprr& alpha, int m )
{
  hprr s2m = sqrt2pow(m);
  return alpha / s2m;
}

halves_t findhalves(const hprr& alpha, int m, const hprr& delta)
{
  hprr s2m = sqrt2pow(m);
  hprr W(weight(alpha,m));
  hprr epsilon = delta * s2m;
  long b = to_long(floor(-s2m));
  long b_max = to_long(ceil(s2m));
  hprr V = alpha * s2m - hprr(b) * hprHelpers::sqrt2();

  halves_t R;

  ft eps(epsilon);
  double w = to_double(W);
  ft v(V);
  ft ms2(-hprHelpers::sqrt2());

  if( to_ld(epsilon) < 0.5 ) // optimized version
  {
    while( b <= b_max )
    {
      long a = to_long(round(v));
      R.push_back(make_pair(w*to_double(v-a),b));
      b++;
      v+=ms2;
    }
  }
  else //generic version
  {
    while( b <= b_max )
    {
      long a_min = to_long(ceil(v-eps));
      long a_max = to_long(floor(v+eps));
      for( long a = a_min; a <= a_max; ++a )
        R.push_back(make_pair(w*to_double(v-a),b));
      b++;
      v+=ms2;
    }

  }

  std::sort(R.begin(),R.end());
  return R;
}
