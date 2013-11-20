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


#include "symbolic_angle.h"

#include <tuple>
using namespace std;

symbolic_angle::operator hprr() const
{
  if( pi )
    return hprr( numerator ) / hprr( denominator ) * hprHelpers::pi() ;
  else
    return hprr( numerator ) / hprr( denominator ) ;
}

std::ostream& operator<<( std::ostream& out , const symbolic_angle& a )
{
  out << a.numerator << "," << a.denominator << "," << a.pi;
  return out;
}

bool symbolic_angle::operator <(const symbolic_angle &a) const
{
  return make_tuple(numerator,denominator,pi) <
      make_tuple(a.numerator,a.denominator,a.pi);
}
