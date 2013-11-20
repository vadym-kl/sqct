//     Copyright (c) 2012 Vadym Kliuchnikov sqct(dot)software(at)gmail(dot)com, Dmitri Maslov, Michele Mosca
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

#ifndef HPRHELPERS_H
#define HPRHELPERS_H

#include <gmpxx.h>

#define MPFR_REAL_DATA_PUBLIC
#include "real.hpp"
#include <complex>

typedef mpfr::real<256> hprr;

/// \brief Helper functions and data for high precision arithmetic
struct hprHelpers
{
    typedef mpfr::real<256> hpr_real;
    /// \brief High precision real type used in the project
    /// \brief High precision complex type used in the project
    typedef std::complex<hpr_real> hpr_complex;
    /// \brief Transforms high precision complex number into machine complex
    static void convert( const hpr_complex& from, std::complex<double>& to );
    /// \brief Transforms high precision real number into double
    static void convert( const hpr_real& from, double& to );
    /// \brief Transforms high precision real number into double
    static double toMachine( const hpr_real& from );
    /// \brief Transforms high precision complex number into machine complex
    static std::complex<double> toMachine( const hpr_complex& from );

    /// \brief High precision \f$ \pi \f$
    static const hpr_real& pi();
    /// \brief High precision one
    static const hpr_real& one();
    /// \brief High precision one
    static const hpr_real& two();
    /// \brief High precision 1/2
    static const hpr_real& half();
    /// \brief High precision -1/2
    static const hpr_real& mhalf();
    /// \brief High precision \f$ \frac{\sqrt{2}}{2} \f$
    static const hpr_real& sqrt2ov2();
    /// \brief High precision \f$ \frac{\sqrt{2}}{2} \f$
    static const hpr_real& sqrt2();
};

hprr pow2( int n );
long double to_ld( const hprr& a );

long to_long( const hprr& a );
mpz_class to_mpz( const hprr& a );
double to_double( const hprr& a );

hprr sqrt2pow( long p );


#endif // HPRHELPERS_H
