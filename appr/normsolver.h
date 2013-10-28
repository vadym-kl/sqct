//     Copyright (c) 2013 Vadym Kliuchnikov sqct(dot)software(at)gmail(dot)com
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
//     Based on "Practical approximation of single-qubit unitaries by single-qubit quantum Clifford and T circuits"
//     by Vadym Kliuchnikov, Dmitri Maslov, Michele Mosca
//     http://arxiv.org/abs/1212.6964v1, new version of the paper describing this modification: [to be published]

#ifndef NORMSOLVER_H
#define NORMSOLVER_H

#include "gatelibrary.h"
#include "rint.h"
#include "factorzs2.h"

#include <pari/pari.h>

/// \brief Solves norm equation |z|^2 = A + \sqrt{2} B, where z is from Z[exp(i pi/4)].
/// PARI (http://en.wikipedia.org/wiki/PARI/GP) C interface is used for this.
class normSolver
{
public:
  typedef matrix2x2<mpz_class> m;
  normSolver();

  /// \brief rsh specifies A + \sqrt{2}B and res contains one of possible solutions to the norm equation
  bool solve( const ring_int_real<mpz_class>& rhs, ring_int<mpz_class>& res ) const;

  bool solve( const mpz_class& rhs, ring_int_real<mpz_class>& res ) const;

  /// \brief Uses norm equation solver to find an entry y of the unitary
  /// 1/2^denompower * {{u00,-y^{\dagger}},{y,u00^{\dagger}}};
  /// \return true if such a unitary exists and writes it down into matr
  bool solve( const ring_int<mpz_class>& u00, int denompower, m& matr ) const;

  zfactorization factor( const mpz_class& number ) const;

  /// \brief Instance of normSolver to be used for all calls; allows to avoid multiple initialization of PARI
  static const normSolver& instance();
private:
  GEN rnf; ///< PARI object for the extension Z[exp(i pi/4)] / Z[sqrt{2}]
  GEN zs2; /// Solving norm equations in Z[sqrt{2}]
};

#endif // NORMSOLVER_H
