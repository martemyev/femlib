#include "solver.h"
#include "matrix.h"
#include "vector.h"

NAMESPACE_FEM_OPEN


Solver::Solver(const Matrix &matrix)
{
  KSPCreate(PETSC_COMM_SELF, &_ksp);
  KSPSetOperators(_ksp, matrix.mat(), matrix.mat(), SAME_PRECONDITIONER);
  KSPSetTolerances(_ksp, 1e-12, 1e-30, 1e+5, 10000);
}



void Solver::solve(const Vector &rhs,
                   Vector &solution)
{
  Vec sol_vec;
  VecDuplicate(rhs.vec(), &sol_vec);
  KSPSolve(_ksp, rhs.vec(), sol_vec);
  solution.init(sol_vec);
}



NAMESPACE_FEM_CLOSE
