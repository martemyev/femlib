#include "vector.h"
#include "auxiliary_functions.h"


NAMESPACE_FEM_OPEN


Vector::Vector()
  : _size(0)
{ }



Vector::Vector(const Vector &vec)
  : _size(vec._size)
{
  if (_size)
    VecDuplicate(vec._data, &_data);
}



Vector& Vector::operator =(const Vector &vec)
{
  _size = vec._size;
  if (_size)
    VecDuplicate(vec._data, &_data);
  return *this;
}



Vector::Vector(int size)
  : _size(size)
{
  VecCreateSeq(PETSC_COMM_SELF, _size, &_data);
}



Vector::~Vector()
{
  VecDestroy(&_data);
}



void Vector::size(int s)
{
  _size = s;
  VecCreateSeq(PETSC_COMM_SELF, _size, &_data);
}



int Vector::size() const
{
  return _size;
}



double Vector::operator ()(int i) const
{
  expect(_size > 0, "Empty vector");
  expect(i >= 0 && i < _size, "Incorrect input parameter");
  double value; // this one will be returned
  // check the returning status ???
  VecGetValues(_data, 1, &i, &value);
  return value;
}



void Vector::values(std::vector<double> &v) const
{
  expect(_size > 0, "Empty vector");
  std::vector<int> idx(_size);
  for (int i = 0; i < _size; ++i)
    idx[i] = i; // idx = { 0, 1, 2, 3, .... }
  v.resize(_size);
  VecGetValues(_data, _size, &idx[0], &v[0]); // extract the values of the vector
}



void Vector::add(int i, double value)
{
  expect(_size > 0, "Empty vector");
  expect(i >= 0 && i < _size, "Incorrect input parameter");
  VecSetValue(_data, i, value, ADD_VALUES);
}



void Vector::insert(int i, double value)
{
  expect(_size > 0, "Empty vector");
  expect(i >= 0 && i < _size, "Incorrect input parameter");
  VecSetValue(_data, i, value, INSERT_VALUES);
}



const Vec& Vector::vec() const
{
  expect(_size > 0, "Empty vector");
  return _data;
}



void Vector::init(const Vec &vector)
{
  VecDuplicate(vector, &_data);
  VecGetSize(_data, &_size);
}

NAMESPACE_FEM_CLOSE
