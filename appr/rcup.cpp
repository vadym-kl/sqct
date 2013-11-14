#include "rcup.h"
#include "findhalves.h"

#include <array>
#include <algorithm>

using namespace std;

double min_positive_first( const halves_t& hv )
{
  halves_t::value_type z(0,0);
  auto r = std::upper_bound(hv.begin(),hv.end(),z);
  if( r != hv.end() )
    return r->first;
  else
    return 1.;
}

zwt to_zwt( const std::array<long,5>& val )
{
  return zwt(val[0],val[1]+val[3],val[2],val[1]-val[3]);
}

typedef std::pair< halves_t::const_iterator, halves_t::const_iterator> it_interval_t;

it_interval_t get_interval(const halves_t& hv, double min_v, double max_v )
{
  it_interval_t r;
  halves_t::value_type l(min_v,0);
  halves_t::value_type h(max_v,0);
  r.first = std::upper_bound(hv.begin(),hv.end(),l);
  if( r.first != hv.begin() ) --r.first;
  r.second = std::upper_bound(hv.begin(),hv.end(),h);
  return r;
}

void rcup::merge_halves( const halves_t& re, const halves_t& im, int k )
{
  size_t sz = re.size();

  for( size_t j = 0; j < sz; ++j )
  {
    const auto& re_j = re[j];
    long a = recover_real(re_j,k);
    auto II = get_interval( im, I.first - re_j.first,I.second - re_j.first );
    for( auto jj = II.first; jj != II.second; ++jj )
    {
      long c = recover_imag(*jj,k);;
      out.push_back(make_pair(jj->first + re_j.first,std::array<long,5>{a,re_j.second,c,jj->second,k}));
    }
  }
}

long rcup::recover_real(const halves_t::value_type &a, int k) const
{
  hprr v = cos2m[k] - a.second * hprHelpers::sqrt2();
  hprr v_minus_a = hprr(a.first) / reW[k];
  return to_long(v - v_minus_a);
}

long rcup::recover_imag(const halves_t::value_type &a, int k) const
{
  hprr v = sin2m[k] - a.second * hprHelpers::sqrt2();
  hprr v_minus_a = hprr(a.first) / imW[k];
  return to_long(v - v_minus_a);
}

rcup::rcup(long n, const hprr &phi, const hprr &delta)
{
  rcup_res R;
  R.first = delta;

  long m = floor(n/2) + 2;
  halves_t L_re[2];
  halves_t L_im[2];
  for( int k = 0; k < 2; ++k )
  {
    cos2m[k] = cos( phi-hprHelpers::pi()*k/hprr(8) ) * sqrt2pow(m);
    sin2m[k] = sin( phi-hprHelpers::pi()*k/hprr(8) ) * sqrt2pow(m);
    reW[k] = weight( cos( phi-hprHelpers::pi()*k/hprr(8)) ,m );
    imW[k] = weight( sin( phi-hprHelpers::pi()*k/hprr(8)) ,m );

    L_re[k] = findhalves(cos(phi-hprHelpers::pi()*k/hprr(8) ), m, delta  );
    L_im[k] = findhalves(sin(phi-hprHelpers::pi()*k/hprr(8) ), m, delta  );
  }

  I = interval_t( 0,min( min_positive_first(L_re[0]), min_positive_first(L_re[1]) ) );
  // Question 1: is doulbe precision sufficient for the purpose of this part ? -- assume yes
  // Question 2: how likely is to get repetitions of epsilon -- assume that it is not likely

  while( I.first < delta )
  {
    res_tuples_arr out;
    for( int k = 0; k < 2; ++k )
      merge_halves( L_re[k], L_im[k], k);

    std::sort(out.begin(),out.end());

    size_t sz_max = out.size();

    for( size_t j = 0; j < sz_max; ++j )
    {
      const auto& o = out[j];
      zwt xp = to_zwt(o.second);

      R.second = min_t_count(xp,m,o.second[4]);
      if( n == R.second.min_t_count )
      {
        R.first = o.first;
        return;
      }
    }

    double t = 2. * I.second - I.first;
    I.first = I.second;
    I.second = t;
  }

  R.second.y.clear();
}

