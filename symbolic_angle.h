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


#ifndef SYMBOLIC_ANGLE_H
#define SYMBOLIC_ANGLE_H

#include "hprhelpers.h"

struct symbolic_angle
{
  long numerator;
  long denominator;
  bool pi;
  operator hprr() const;
  bool operator < ( const symbolic_angle& angle ) const;
};

std::ostream& operator<<( std::ostream& out , const symbolic_angle& a );

#endif // SYMBOLIC_ANGLE_H
