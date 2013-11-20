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


#include "fixedpoint.h"

#include <sstream>

using namespace std;
using namespace boost::multiprecision;


fixedpoint::fixedpoint(const hprr &val)
{
  stringstream ss;
  ss << to_mpz(pow2(frac_bits) * val);
  m_value = bigint_t(ss.str());
}

fixedpoint::fixedpoint(const fixedpoint &val) : m_value(val.m_value)
{
}

fixedpoint fixedpoint::operator+(const fixedpoint &rhs) const
{
  fixedpoint r(*this);
  r.m_value+= rhs.m_value;
  return r;
}

fixedpoint fixedpoint::operator-(const fixedpoint &rhs) const
{
  fixedpoint r(*this);
  r.m_value-= rhs.m_value;
  return r;
}

fixedpoint fixedpoint::operator-(long rhs) const
{
  fixedpoint r(*this);
  bigint_t rh(rhs);
  rh <<= frac_bits;
  r.m_value-= rhs;
  return r;
}

fixedpoint &fixedpoint::operator+=(const fixedpoint &rhs)
{
  m_value+= rhs.m_value;
  return *this;
}

double frac_to_double( const fixedpoint::bigint_t& frac )
{
  fixedpoint::bigint_t fr(frac);
  fr >>= 64;
  long high = fr.convert_to<long>();
  fr <<= 64;
  long low = (frac-fr).convert_to<long>();
  long double r = ldexp((long double)(low),-128);
  if( high > 0 )
  {
    r += ldexp((long double)(high),-64);
  }
  return r;
}

std::pair<long, double> fixedpoint::round_ex() const
{
  std::pair<long, double> res;
  bigint_t int_part(abs(m_value));
  bigint_t frac_part(int_part);

  int_part >>= frac_bits;
  res.first = int_part.convert_to<long>();
  int_part <<= frac_bits;
  frac_part -= int_part;

  bigint_t t(frac_part);
  t >>= (frac_bits - 1);

  bool neg = false;
  if( t > 0 )
  {
    res.first+=1;
    frac_part-=one;
    neg = true;
  }

  res.second = frac_to_double(abs(frac_part));
  if( neg ) res.second = -res.second;

  if( m_value < 0 )
  {
    res.first = -res.first;
    res.second = -res.second;
  }

  return res;
}


double fixedpoint::to_double() const
{
  auto res = round_ex();
  return double(res.first) + res.second;
}

long ceil(const fixedpoint &val)
{
  auto res = val.round_ex();
  if( res.second <= 0. )
    return res.first;
  else
    return res.first + 1;
}


long floor(const fixedpoint &val)
{
  auto res = val.round_ex();
  if( res.second < 0. )
    return res.first - 1;
  else
    return res.first;
}

double to_double(const fixedpoint &val)
{
  return val.to_double();
}

//fixedpoint::bigint_t fixedpoint::mask("0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF");
//fixedpoint::bigint_t fixedpoint::half("0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF");
fixedpoint::bigint_t fixedpoint::one("0x100000000000000000000000000000000");

