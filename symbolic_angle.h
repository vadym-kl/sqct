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
