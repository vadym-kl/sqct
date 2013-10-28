#include "test.h"

#include "appr/normsolver.h"
#include "solvenormequation.h"
#include "output.h"

#include <iostream>
#include <sstream>
#include <cassert>

using namespace std;

template < class T >
bool test_string( const T& val, const string& str )
{
  stringstream ss; ss << val;
  bool eq = (ss.str() == str);
  if( !eq )
  {
    cout << "got:" << ss.str() << endl;
    cout << "expected:" << str << endl;
  }
  return eq;
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

void zs2normEquationTest()
{
  {
    mpz_class v = 7;
    ring_int_real<mpz_class> ans;
    bool r = normSolver::instance().solve(v,ans);
    assert(r);
    assert(test_string(ans,
                       "zs2[3,-1]"));
  }

  {
    mpz_class v = 11;
    ring_int_real<mpz_class> ans;
    bool r = normSolver::instance().solve(v,ans);
    assert(!r);
  }
}

void unit_log_test()
{
  {
    zs2type u00(-41,29);
    auto r = unit_log(u00);
    assert( unit_power<mpz_class>(r) == u00 );
    assert( test_string(r,"{1,5}") );
  }

  {
    zs2type u00(41,-29);
    auto r = unit_log(u00);
    assert( unit_power<mpz_class>(r) == u00 );
    assert( test_string(r,"{-1,5}") );
  }

  {
    zs2type u00(99,70);
    auto r = unit_log(u00);
    assert( unit_power<mpz_class>(r) == u00 );
    assert( test_string(r,"{1,-6}") );
  }

  {
    zs2type u00(-99,-70);
    auto r = unit_log(u00);
    assert( unit_power<mpz_class>(r) == u00 );
    assert( test_string(r,"{-1,-6}") );
  }
}

void zs2FactoringTest()
{
  {
    ring_int_real<mpz_class> num7(91,11);
    auto n = num7*num7*(num7.g_conjugate());
    auto res = factorize( n );
    assert(zs2type(res) == n );
    assert(test_string(res.prime_factors,
                       "{{zs2[91,11],2},{zs2[91,-11],1}}"));
  }

  {
    ring_int_real<mpz_class> num1(109,40);
    auto n = num1*num1*(num1.g_conjugate());
    auto res = factorize( n );
    assert(zs2type(res) == n );
    assert(test_string(res.prime_factors,
                       "{{zs2[109,-40],1},{zs2[109,40],2}}"));
  }

  {
    ring_int_real<mpz_class> num1(-109,40);
    auto n = num1*num1*(num1.g_conjugate());
    auto res = factorize( n );
    assert(zs2type(res) == n );
  }
}

bool norm_solver_sub_test( const zs2type& rhs )
{
  auto r = solve_norm_equation(rhs);
  if( r.exists )
  {
    bool eq = zwt(r).abs2() == rhs;
    if( !eq )
    {
      cout << "given: " << rhs << endl;
      cout << "got  : " << zwt(r).abs2() << endl;
    }
    assert( eq );
  }
  return r.exists;
}

void norm_solver_test()
{
  zs2type ram(2,1);
  zs2type num1(109,40);
  zs2type num7(91,11);
  zs2type num3(8011,0);
  zs2type num5(8053,0);
  assert(norm_solver_sub_test(ram));
  assert(norm_solver_sub_test(num1));
  assert(!norm_solver_sub_test(num7));
  assert(!norm_solver_sub_test(num7*ram));
  assert(norm_solver_sub_test(num7*num7*ram));
  assert(norm_solver_sub_test(num3));
  assert(norm_solver_sub_test(num5));
  assert(norm_solver_sub_test(num7*num7*ram*ram*ram*num1*num3*num5));
}

void run_tests()
{
  cout << "Testing started" << endl;
  z_factoring_test();
  zs2normEquationTest();
  zs2FactoringTest();
  unit_log_test();
  norm_solver_test();
  cout << "Testing finished" << endl;
}
