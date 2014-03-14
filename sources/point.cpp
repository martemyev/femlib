#include "point.h"
#include "auxiliary_functions.h"
#include <vector>

NAMESPACE_FEM_OPEN



Point::Point()
{
  for (unsigned i = 0; i < n_coord; ++i)
    _coord[i] = 0.;
}



Point::Point(const double coordinates[])
{
  for (unsigned i = 0; i < n_coord; ++i)
    _coord[i] = coordinates[i];
}



Point::Point(const double x_coord,
             const double y_coord,
             const double z_coord)
{
  _coord[0] = x_coord;
  if (n_coord > 1) _coord[1] = y_coord;
  if (n_coord > 2) _coord[2] = z_coord;
}



Point::Point(const Point &p)
{
  for (unsigned i = 0; i < n_coord; ++i)
    _coord[i] = p._coord[i];
}



Point& Point::operator =(const Point &p)
{
  for (unsigned i = 0; i < n_coord; ++i)
    _coord[i] = p._coord[i];
  return *this;
}



double Point::coord(unsigned int number) const
{
  expect(number < n_coord,
         "The number of coordinate is incorrect: " +
         d2s(number) + ". It should be in the range: [0, " +
         d2s(n_coord) + ")");

  return _coord[number];
}



void Point::coord(unsigned int number, double value)
{
  expect(number < n_coord,
         "The number of coordinate is incorrect: " +
         d2s(number) + ". It should be in the range: [0, " +
         d2s(n_coord) + ")");

  _coord[number] = value;
}



Point& Point::operator /=(double d)
{
  for (unsigned i = 0; i < n_coord; ++i)
    _coord[i] /= d;
  return *this;
}



Point operator -(const Point &p1, const Point &p2)
{
  Point res;
  for (unsigned i = 0; i < Point::n_coord; ++i)
    res._coord[i] = p1._coord[i] - p2._coord[i];
  return res;
}



std::ostream& operator <<(std::ostream &os, const Point &p)
{
  os << "(";
  for (unsigned i = 0; i < Point::n_coord; ++i)
    os << p._coord[i] << ",";
  os << ")";
  return os;
}


NAMESPACE_FEM_CLOSE
