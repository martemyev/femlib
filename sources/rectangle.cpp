#include "rectangle.h"
#include "auxiliary_functions.h"
#include <cmath>


Rectangle::Rectangle()
{ }



Rectangle::~Rectangle()
{ }



// A rectangle has the following numeration of vertices
// 2 ----- 3
// |       |
// |       |
// 0 ----- 1
// Therefore the x-coordinates of the (0-th and 2-nd) and (1-st and 3-rd) vertices
// and y-coordinates of the (0-th and 1-st) and (2-nd and 3-rd) vertices
// should coincide
Rectangle::Rectangle(const std::vector<unsigned int> &ver,
                     const std::vector<Point> &mesh_vertices,
                     const unsigned int mat_id,
                     const unsigned int part_id,
                     const std::vector<unsigned int> &g_cells)
  : Quadrangle(ver, mesh_vertices, mat_id, part_id, g_cells)
{
  check_rectangle();
}



Rectangle::Rectangle(const unsigned int v1,
                     const unsigned int v2,
                     const unsigned int v3,
                     const unsigned int v4,
                     const std::vector<Point> &mesh_vertices,
                     const unsigned int mat_id,
                     const unsigned int part_id,
                     const std::vector<unsigned int> &g_cells)
  : Quadrangle(v1, v2, v3, v4, mesh_vertices, mat_id, part_id, g_cells)
{
  check_rectangle();
}



Rectangle::Rectangle(const Rectangle &rect)
  : Quadrangle(rect),
    _dofs(rect._dofs)
{ }



Rectangle& Rectangle::operator =(const Rectangle &rect)
{
  Quadrangle::operator =(rect);
  _dofs = rect._dofs;
  return *this;
}



void Rectangle::local_stiffness_matrix(const double coef_beta, double **loc_mat) const
{
  switch (_dofs.size())
  {
  case n_dofs_first:
    {
      const double hx = _X[1] - _X[0]; // the x-length of the cell
      const double hy = _Y[2] - _Y[0]; // the y-length of the cell
      const double m11 = hy/hx + hx/hy;
      const double m12 = -hy/hx + hx/(2.*hy);
      const double m13 = hy/(2.*hx) - hx/hy;
      const double m14 = -hy/(2.*hx) - hx/(2.*hy);
      const double mat[][n_dofs_first] = { { m11, m12, m13, m14 },
                                           { m12, m11, m14, m13 },
                                           { m13, m14, m11, m12 },
                                           { m14, m13, m12, m11 }
                                         };
      for (int i = 0; i < n_dofs_first; ++i)
        for (int j = 0; j < n_dofs_first; ++j)
          loc_mat[i][j] = coef_beta * mat[i][j] / 3.;
      return;
    }
  default:
    require(false, "Unknown order of bf");
  }
}



void Rectangle::local_mass_matrix(const double coef_alpha, double **loc_mat) const
{
  switch (_dofs.size())
  {
  case n_dofs_first:
    {
      const double hx = _X[1] - _X[0]; // the x-length of the cell
      const double hy = _Y[2] - _Y[0]; // the y-length of the cell
      const double mat[][n_dofs_first] = { { 4, 2, 2, 1 },
                                           { 2, 4, 1, 2 },
                                           { 2, 1, 4, 2 },
                                           { 1, 2, 2, 4 }
                                         };
      for (int i = 0; i < n_dofs_first; ++i)
        for (int j = 0; j < n_dofs_first; ++j)
          loc_mat[i][j] = coef_alpha * hx * hy * mat[i][j] / 36.;
      return;
    }
  default:
    require(false, "Unknown order of bf");
  }
}



void Rectangle::local_rhs_vector(double (*rhs_func)(const Point& p, double t),
                                 const std::vector<Point> &points,
                                 const double time,
                                 double *loc_vec) const
{
  switch (_dofs.size())
  {
  case n_dofs_first:
    {
      const double hx = _X[1] - _X[0]; // the x-length of the cell
      const double hy = _Y[2] - _Y[0]; // the y-length of the cell
      const double mat[][n_dofs_first] = { { 4, 2, 2, 1 },
                                           { 2, 4, 1, 2 },
                                           { 2, 1, 4, 2 },
                                           { 1, 2, 2, 4 }
                                         };
      for (int i = 0; i < n_dofs_first; ++i)
      {
        loc_vec[i] = 0;
        for (int j = 0; j < n_dofs_first; ++j)
          loc_vec[i] += hx * hy * mat[i][j] * rhs_func(points[_dofs[j]], time) / 36.;
      }
      return;
    }
  default:
    require(false, "Unknown order of bf");
  }
}



void Rectangle::check_rectangle() const
{
  expect(fabs(_X[0] - _X[2]) < FLOAT_NUMBERS_EQUALITY_TOLERANCE &&
         fabs(_X[1] - _X[3]) < FLOAT_NUMBERS_EQUALITY_TOLERANCE &&
         fabs(_Y[0] - _Y[1]) < FLOAT_NUMBERS_EQUALITY_TOLERANCE &&
         fabs(_Y[2] - _Y[3]) < FLOAT_NUMBERS_EQUALITY_TOLERANCE,
         "This is not a rectangle, since the angles are not right, OR the numeration of the vertices is not the one that was expected."
         " The differences are: fabs(_X[0] - _X[2]) = " + d2s(fabs(_X[0] - _X[2])) +
         " fabs(_X[1] - _X[3]) = " + d2s(fabs(_X[1] - _X[3])) +
         " fabs(_Y[0] - _Y[1]) = " + d2s(fabs(_Y[0] - _Y[1])) +
         " fabs(_Y[2] - _Y[3]) = " + d2s(fabs(_Y[2] - _Y[3])));
}



unsigned int Rectangle::n_dofs() const
{
  return _dofs.size();
}



void Rectangle::n_dofs(unsigned int number)
{
  _dofs.resize(number, 0);
}


void Rectangle::dof(unsigned int number, unsigned int value)
{
  expect(number < _dofs.size(), "");
  _dofs[number] = value;
}



unsigned int Rectangle::dof(unsigned int number) const
{
  expect(number < _dofs.size(), "");
  return _dofs[number];
}
