#include "factorzs2.h"
#include "appr/normsolver.h"

zfactorization factorize(const ztype &val)
{
  return normSolver::instance().factor(val);
}


zs2factorization factorize(const zs2type &val, const zfactorization &factors)
{

}
