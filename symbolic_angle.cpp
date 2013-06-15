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
