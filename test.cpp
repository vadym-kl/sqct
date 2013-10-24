#include "test.h"

#include "appr/normsolver.h"
#include "output.h"

#include <iostream>
#include <sstream>
#include <cassert>

using namespace std;

template < class T >
bool test_string( const T& val, const string& str )
{
  stringstream ss; ss << val;
  return ss.str() == str;
}

void z_factoring_test()
{
  {
    mpz_class v("-4396707932524505161140129847885469217");
    auto r = normSolver::instance().factor(v);
    assert(test_string(r.prime_factors,
                       "{{3,1},{251,1},{3216141233,1},{1815505332050029736149633,1}}"));
  }

  {
    mpz_class v = 32 * 5;
    auto r = normSolver::instance().factor(v);
    assert(test_string(r.prime_factors,
                       "{{2,5},{5,1}}"));
  }
}

void run_tests()
{
  cout << "Testing started" << endl;
  z_factoring_test();
  cout << "Testing finished" << endl;
}
