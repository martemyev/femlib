#ifndef FEM_MATRIX_H
#define FEM_MATRIX_H

#include "config.h"
#include "petscmat.h"

NAMESPACE_FEM_OPEN



/**
 * Wrapper for PETSc Mat class
 */
class Matrix
{
public:
  Matrix();
  Matrix(const Matrix &mat);
  Matrix& operator =(const Matrix &mat);

  virtual ~Matrix();

            /**
             * Get the i-th element
             * @param i - the number of the row of the element we want to return
             * @param j - the number of the column of the element we want to return
             */
  virtual double operator ()(unsigned i, unsigned j) = 0;

  const Mat& mat() const;

private:
  unsigned int _n_rows;
  unsigned int _n_cols;
  Mat _data;
};


NAMESPACE_FEM_CLOSE

#endif // FEM_MATRIX_H
