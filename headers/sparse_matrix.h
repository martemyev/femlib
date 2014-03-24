#ifndef FEM_SPARSE_MATRIX_H
#define FEM_SPARSE_MATRIX_H

#include "config.h"
#include "matrix.h"

NAMESPACE_FEM_OPEN

class Pattern;


class SparseMatrix
{
public:
  SparseMatrix(const Pattern &pattern);

private:


};


NAMESPACE_FEM_CLOSE


#endif // FEM_SPARSE_MATRIX_H
