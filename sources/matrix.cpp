#include "matrix.h"
#include "auxiliary_functions.h"

NAMESPACE_FEM_OPEN


Matrix::Matrix()
  : _n_rows(0),
    _n_cols(0)
{ }



Matrix::Matrix(const Matrix &mat)
  : _n_rows(mat._n_rows),
    _n_cols(mat._n_cols)
{
  if (_n_rows && _n_cols)
    MatDuplicate(mat._data, MAT_COPY_VALUES, &_data);
}



Matrix& Matrix::operator =(const Matrix &mat)
{
  _n_rows = mat._n_rows;
  _n_cols = mat._n_cols;
  if (_n_rows && _n_cols)
    MatDuplicate(mat._data, MAT_COPY_VALUES, &_data);
  return *this;
}



Matrix::~Matrix()
{
  MatDestroy(&_data);
}



const Mat& Matrix::mat() const
{
  expect(_n_rows && _n_cols, "The matrix is empty");
  return _data;
}


NAMESPACE_FEM_CLOSE
