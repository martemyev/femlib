#ifndef FEM_SOLVER_H
#define FEM_SOLVER_H

#include "config.h"
#include "petscksp.h"

NAMESPACE_FEM_OPEN

class Matrix;
class Vector;



class Solver
{
public:
  Solver(const Matrix &matrix);

  void solve(const Vector &rhs,
             Vector &solution);

private:
  KSP _ksp;

};

NAMESPACE_FEM_CLOSE


#endif // FEM_SOLVER_H
