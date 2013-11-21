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

#include "normsolver.h"
#include <iostream>

#include <cassert>

#define PARI_DEBUG

using namespace std;

normSolver::normSolver()
{
  pari_init_opts(40000000l,1l << 24,INIT_DFTm);
  rnf = gp_read_str("rnfisnorminit(z^2-2,x^2+1)");
  zs2 = gp_read_str("bnfinit(z^2-2)");
}

static void genToMpz( GEN gen, mpz_class& out )
{
  if( gen != 0 )
  {
    char* gen_str = GENtostr(gen);
    out = mpz_class(gen_str);
    pari_free(gen_str);
  }
  else
    out = 0;
}

static void genToMpz( GEN gen, mpz_class& out, bool& div )
{
  div = false;
  if( gen != 0 )
  {
    if( typ(gen) == t_INT )
    {
      char* gen_str = GENtostr(gen);
      out = mpz_class(gen_str);
      pari_free(gen_str);
    }
    else if( typ(gen) == t_FRAC )
    {
      char* gen_str = GENtostr(gel(gen,1));
      out = mpz_class(gen_str);
      pari_free(gen_str);
      div = true;
      if( ! gequalgs(gel(gen,2),2) )
        throw std::logic_error("unexpected denominator");
    }
    else
      throw std::logic_error("wrong data type");
  }
  else
    out = 0;
}

void gen_to_pair( GEN re_part, GEN& re_a, GEN& re_b )
{
  if( typ(re_part) == t_POLMOD )
  {

#ifdef PARI_DEBUG
auto lgsln2 = lg(re_part);
if( lgsln2 != 3 )
  throw std::logic_error(__FILE__":6");
#endif

    if( typ(gel(re_part,2)) == t_POL )
    {
      GEN re = gel(re_part,2);

#ifdef PARI_DEBUG
auto lgre = lg(re_part);
if( lgre != 3 && lgre != 4 )
  throw std::logic_error(__FILE__":6");
#endif
      if( lg(re) >= 3 )
        re_a = gel(re,2);
      if( lg(re) == 4 )
        re_b = gel(re,3);
    }
    else if( typ(gel(re_part,2)) == t_INT )
    {
      re_a = gel(re_part,2);
    }
    else
      throw std::logic_error(__FILE__":re2:unexpected");
  }
  else if( typ(re_part) == t_INT )
  {
    re_a = re_part;
  }
  else
    throw std::logic_error(__FILE__":re:unexpected");
}

bool normSolver::solve(const ring_int_real<mpz_class>& rhs, ring_int<mpz_class> &res) const
{

  pari_sp av = avma;

  stringstream ss;
  ss << "(" << rhs[0] << ")+(" << rhs[1] << ")*z";
  GEN in = gp_read_str(ss.str().c_str());
  GEN solution = rnfisnorm(rnf,in,0);

#ifdef PARI_DEBUG
  if( lg(solution) != 3 )
    throw std::logic_error(__FILE__":1");
#endif

  bool ok =  gequal1(gel(solution,2));
  if( !ok ) //there is no solution
  {
    avma = av;
    return false;
  }

#ifdef PARI_DEBUG
  if( typ(gel(solution,1)) != t_POLMOD )
    throw std::logic_error(__FILE__":2");
  if( lg(gel(solution,1)) != 3 )
    throw std::logic_error(__FILE__":3");
#endif

  GEN re_a = 0, re_b = 0, im_a = 0, im_b = 0;
  GEN sln = gel(gel(solution,1),2);



#ifdef PARI_DEBUG
  if( typ(sln) != t_POL )
    throw std::logic_error(__FILE__":4");
  auto lgsln = lg(sln);
  if( lgsln != 3 && lgsln != 4 )
    throw std::logic_error(__FILE__":5");
#endif


  if( lg(sln) >= 3 )
    gen_to_pair( gel(sln,2), re_a, re_b );

  if( lg(sln) == 4 )
    gen_to_pair( gel(sln,3), im_a, im_b );

  genToMpz( re_a, res[0] );
  genToMpz( im_a, res[2] );

  bool re_b_fr, im_b_fr;
  mpz_class re_bz, im_bz;
  genToMpz( re_b, re_bz, re_b_fr );
  genToMpz( im_b, im_bz, im_b_fr );

  if( re_b_fr && im_b_fr )
  {
    res[1] = re_bz + im_bz;
    res[1] /= 2;
    res[3] = im_bz - re_bz;
    res[3] /= 2;
  }
  else if( !im_b_fr && !re_b_fr )
  {
    res[1] = re_bz + im_bz;
    res[3] = im_bz - re_bz;
  }
  else
  {
    cout << ss.str().c_str() << endl;
    pari_printf("(%Ps)+Sqrt[2]*(%Ps)+I*((%Ps)+Sqrt[2]*(%Ps))\n",re_a,re_b,im_a,im_b);
    throw std::logic_error("unexpected solution");
  }

  avma = av;
  return true;
}

bool normSolver::solve(const mpz_class &rhs, ring_int_real<mpz_class> &res) const
{
   pari_sp av = avma;

   stringstream ss;
   ss << rhs;
   GEN in = gp_read_str(ss.str().c_str());
   GEN sln = bnfisnorm(zs2,in,1);

   assert( lg(sln) == 3 );

   bool ok =  gequal1(gel(sln,2));
   if( !ok ) //there is no solution
   {
     avma = av;
     return false;
   }

   GEN val = gel(sln,1);
   assert( typ(val) == t_POLMOD );
   GEN val_basis = algtobasis(zs2,val);
   assert( typ(val_basis) == t_COL );
   assert( lg(val_basis) == 3 );
   assert( typ(gel(val_basis,1)) == t_INT );
   assert( typ(gel(val_basis,2)) == t_INT );

   genToMpz( gel(val_basis,1), res[0] );
   genToMpz( gel(val_basis,2), res[1] );

   res[2] = 0;
   res[3] = -res[1];

   res.make_positive();

   avma = av;
   return true;
}

bool normSolver::solve( const ring_int<mpz_class>& u00, int denompower, m& matr ) const
{
  mpz_class pow2 = 1;
  pow2 <<= denompower * 2;
  ring_int_real<mpz_class> norm(u00.abs2());
  if( norm[0] > pow2 )
    return false;
  norm[0] = pow2 - norm[0];
  norm[1] = -norm[1];
  norm[3] = -norm[3];

  ring_int<mpz_class> u10;
  if( !solve(norm,u10) )
    return false;

  matr.set(u00,-u10.conjugate(),
        u10, u00.conjugate(), denompower * 2);

  matr.reduce();
  return true;
}

zfactorization normSolver::factor(const mpz_class &number) const
{
  mpz_class n;
  zfactorization res;

  pari_sp av = avma;

  if( number < 0 )
  {
    n = -number;
    res.sign = -1;
  }
  else
  {
    n = number;
    res.sign = 1;
  }

  GEN in = gp_read_str(n.get_str().c_str());
  assert( typ(in) == t_INT );
  GEN F = Z_factor(in);
  assert( typ(F) == t_MAT );
  assert( lg(F) == 3 );

  GEN primes = gel(F,1);
  assert( typ(primes) == t_COL );

  GEN powers = gel(F,2);
  assert( typ(powers) == t_COL );

  int factors = lg(primes) - 1;
  assert( lg(powers) == factors + 1 );

  for( int i = 0; i < factors; ++i )
  {
    GEN prime = gel(primes,i+1);
    GEN power = gel(powers,i+1);
    assert(typ(prime) == t_INT);
    assert(typ(power) == t_INT);
    mpz_class pr,pw;
    genToMpz(prime,pr);
    genToMpz(power,pw);
    res.prime_factors.push_back(make_pair(pr,pw.get_si()));
  }

  avma = av;

  return res;
}

const normSolver &normSolver::instance()
{
  static normSolver inst;
  return inst;
}

