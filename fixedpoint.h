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


#ifndef FIXEDPOINT_H
#define FIXEDPOINT_H

#include "hprhelpers.h"

#include <ttmath/ttmath.h>

#include <cmath>

#include <boost/multiprecision/cpp_int.hpp>

using boost::multiprecision::backends::cpp_int_backend;
using boost::multiprecision::signed_magnitude;
using boost::multiprecision::unchecked;



struct fixedpoint
{
  fixedpoint( const hprr& val );
  fixedpoint( const fixedpoint& val );

  fixedpoint operator+ ( const fixedpoint& rhs ) const;
  fixedpoint operator- ( const fixedpoint& rhs ) const;
  fixedpoint operator- ( long rhs ) const;
  fixedpoint& operator+=( const fixedpoint& rhs );

  std::pair<long,double> round_ex() const;

  double to_double() const;

  static const int bitsize = 192;
  static const int int_bits = 64;
  static const int frac_bits = bitsize - int_bits;

  typedef boost::multiprecision::number< cpp_int_backend<bitsize, bitsize, signed_magnitude, unchecked, void> > bigint_t;

  static bigint_t one;

  bigint_t m_value;

};



long ceil( const fixedpoint& val );

long floor( const fixedpoint& val );

double to_double( const fixedpoint& val );

struct grid_iterator
{
  static const int frac_words = 2;
  static const int int_words = 1;

  typedef ttmath::UInt<frac_words+int_words> uint_type;


  uint_type m_current; //  b*step_offset + 2^{frac_words*64}*m_big_offset

  long m_big_offset;
  uint_type m_offset;

  uint_type m_down_threshold;
  uint_type m_up_threshold;

  grid_iterator( const hprr& initial, const hprr& step_offset, const hprr& threshold, long m )
  {
    int rescale_pow = frac_words*sizeof(long)*8;
    assert( (m/2+ 10) < 63 );
    m_big_offset = 1L << (m/2+10);
    mpz_class init = to_mpz( ldexp( hprr(m_big_offset) + initial, rescale_pow ) );
    //assert( init >= 0 );
    if( init < 0 )
    {
      std::cout << init << std::endl;
      throw std::logic_error("fail");
    }
    m_current = uint_type( init.get_str() );
    m_offset = uint_type( to_mpz( ldexp( step_offset, rescale_pow) ).get_str());
    m_down_threshold = uint_type( to_mpz( ldexp( threshold, rescale_pow) ).get_str() );
    uint_type one(uint_type(1) << rescale_pow);
    m_up_threshold = one - m_down_threshold;

    if( m_up_threshold == one )
    {
      m_up_threshold = one - uint_type(1);
    }

    assert(m_down_threshold >= 0 && m_down_threshold < one );
    assert(m_up_threshold >= 0 && m_up_threshold < one );
  }

  void operator++()
  {
    m_current -= m_offset;
  }

  bool up_close() const
  {
    for( int i = frac_words - 1; i >= 0; --i )
      if( m_current.table[i] < m_up_threshold.table[i] )
        return false;
      else if( m_current.table[i] > m_up_threshold.table[i] )
        return true;
    return true;
  }

  bool down_close() const
  {
    for( int i = frac_words - 1; i >= 0; --i )
      if( m_current.table[i] > m_down_threshold.table[i] )
        return false;
      else if( m_current.table[i] < m_up_threshold.table[i] )
        return true;
    return true;
  }

  long up_a()
  {
    return m_current.table[frac_words] + 1 - m_big_offset;
  }

  long down_a()
  {
    return m_current.table[frac_words] - m_big_offset;
  }

  long double up_frac_part()
  {
    uint_type one(uint_type(1) << (frac_words*sizeof(long)*8) );
    uint_type fr(0);

    for( int i = 0; i < frac_words; ++i )
      fr.table[i] = m_current.table[i];

    one -= fr;

    long double res = 0.;
    for( int i = 0; i < frac_words; ++i )
    {
      if( m_current.table[i] != 0 )
        res += ldexp((long double)(one.table[i]), -(frac_words-i)*sizeof(long)*8 );
    }
    assert( ! isnan(res) );
    return -res;
  }

  long double down_frac_part()
  {
    long double res = 0.;
    for( int i = 0; i < frac_words; ++i )
    {
      if( m_current.table[i] != 0 )
        res += ldexp((long double)(m_current.table[i]), -(frac_words-i)*sizeof(long)*8 );
    }
    assert( ! isnan(res) );
    return res;
  }
};

#endif // FIXEDPOINT_H
