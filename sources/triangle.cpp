#include "triangle.h"
#include "auxiliary_functions.h"
#include "edge.h"
#include "dof_handler.h"
#include "math_functions.h"


Triangle::Triangle()
  : MeshElement(n_vertices, n_edges, n_faces, gmsh_el_type),
    _detD(0.)
{
  for (int i = 0; i < n_vertices; ++i)
    _A[i] = _B[i] = _C[i] = 0.;
  _dis_edges.resize(n_edges, 0);
}



Triangle::~Triangle()
{
  _ghost_cells.clear();
  _dofs.clear();
  _dis_edges.clear();
}



Triangle::Triangle(const std::vector<unsigned int> &ver,
                   const std::vector<Point> &mesh_vertices,
                   const unsigned int mat_id,
                   const unsigned int part_id,
                   const std::vector<unsigned int> &g_cells)
  : MeshElement(n_vertices, n_edges, n_faces, gmsh_el_type)
{
  expect(_vertices.size() == ver.size(),
         "The size of list of vertices (" + d2s(ver.size()) +
         ") is not equal to really needed number of vertices (" + d2s(n_vertices) + ")");
  _vertices = ver;
  _material_id = mat_id;
  _partition_id = part_id;
  _ghost_cells = g_cells;
  _dis_edges.resize(n_edges, 0);

  if (!mesh_vertices.empty())
    calculate_barycentric(mesh_vertices);
}



Triangle::Triangle(const unsigned int v1,
                   const unsigned int v2,
                   const unsigned int v3,
                   const std::vector<Point> &mesh_vertices,
                   const unsigned int mat_id,
                   const unsigned int part_id,
                   const std::vector<unsigned int> &g_cells)
  : MeshElement(n_vertices, n_edges, n_faces, gmsh_el_type)
{
  _vertices[0] = v1;
  _vertices[1] = v2;
  _vertices[2] = v3;
  _material_id = mat_id;
  _partition_id = part_id;
  _ghost_cells = g_cells;
  _dis_edges.resize(n_edges, 0);

  if (!mesh_vertices.empty())
    calculate_barycentric(mesh_vertices);
}



Triangle::Triangle(const Triangle &tri)
  : MeshElement(tri),
    _ghost_cells(tri._ghost_cells),
    _dofs(tri._dofs),
    _detD(tri._detD)
{
  for (int i = 0; i < n_vertices; ++i)
  {
    _A[i] = tri._A[i];
    _B[i] = tri._B[i];
    _C[i] = tri._C[i];
  }
  _dis_edges = tri._dis_edges;
}



Triangle& Triangle::operator =(const Triangle &tri)
{
  MeshElement::operator =(tri);
  _ghost_cells = tri._ghost_cells;
  _dofs = tri._dofs;
  _detD = tri._detD;
  for (int i = 0; i < n_vertices; ++i)
  {
    _A[i] = tri._A[i];
    _B[i] = tri._B[i];
    _C[i] = tri._C[i];
  }
  _dis_edges = tri._dis_edges;
}



unsigned int Triangle::n_ghost_cells() const
{
  return _ghost_cells.size();
}


unsigned int Triangle::ghost_cell(unsigned int num) const
{
  expect(num >= 0 && num < _ghost_cells.size(), "Incorrect input parameter");
  return _ghost_cells[num];
}



unsigned int Triangle::n_dofs() const
{
  return _dofs.size();
}



unsigned int Triangle::dof(unsigned int num) const
{
  expect(num >= 0 && num < _dofs.size(), "Incorrect input parameter");
  return _dofs[num];
}



void Triangle::n_dofs(unsigned int number)
{
  _dofs.resize(number);
}



void Triangle::dof(unsigned int num, unsigned int value)
{
  expect(num >= 0 && num < _dofs.size(), "Incorrect input parameter");
  _dofs[num] = value;
}



void Triangle::local_mass_matrix(const double coef_alpha, double **loc_mat) const
{
  expect(fabs(_detD) > 1e-15,
         "An attempt to calculate a local matrix for a singular triangle (detD = " + d2s(_detD, true) + ")");

  switch (_dofs.size())
  {
  case n_dofs_first:
  {
    const double mat[][n_dofs_first] = { { 2., 1., 1. },
                                         { 1., 2., 1. },
                                         { 1., 1., 2. }
                                       };
    for (int i = 0; i < n_dofs_first; ++i)
      for (int j = 0; j < n_dofs_first; ++j)
        loc_mat[i][j] = coef_alpha * fabs(_detD) * mat[i][j] / 24.;
    return;
  }
  default:
    require(false, "Unknown fe order for generating local matrix");
  }
}



void Triangle::local_stiffness_matrix(const double coef_beta, double **loc_mat) const
{
  expect(fabs(_detD) > 1e-15,
         "An attempt to calculate a local matrix for a singular triangle (detD = " + d2s(_detD, true) + ")");

  switch (_dofs.size())
  {
  case n_dofs_first:
  {
    for (int i = 0; i < n_dofs_first; ++i)
      for (int j = 0; j < n_dofs_first; ++j)
        loc_mat[i][j] = coef_beta * fabs(_detD) * (_A[i]*_A[j] + _B[i]*_B[j]) / 2.;
    return;
  }
  default:
    require(false, "Unknown fe order for generating local matrix");
  }
}



void Triangle::calculate_barycentric(const std::vector<Point> &mesh_vertices)
{
  double x[n_vertices], y[n_vertices]; // coordinates of the triangle's vertices
  for (int i = 0; i < n_vertices; ++i)
  {
    const unsigned int vertex = _vertices[i];
    x[i] = mesh_vertices[vertex].coord(0); // x-coordinate
    y[i] = mesh_vertices[vertex].coord(1); // y-coordinate
  }
  _detD = (x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0]);
  _A[0] = (y[1] - y[2]) / _detD;
  _A[1] = (y[2] - y[0]) / _detD;
  _A[2] = (y[0] - y[1]) / _detD;
  _B[0] = (x[2] - x[1]) / _detD;
  _B[1] = (x[0] - x[2]) / _detD;
  _B[2] = (x[1] - x[0]) / _detD;
  _C[0] = (x[1]*y[2] - x[2]*y[1]) / _detD;
  _C[1] = (x[2]*y[0] - x[0]*y[2]) / _detD;
  _C[2] = (x[0]*y[1] - x[1]*y[0]) / _detD;
}



void Triangle::local_rhs_vector(double(*rhs_func)(const Point &point, double t),
                                const std::vector<Point> &points,
                                const double time,
                                double *loc_vec) const
{
  expect(fabs(_detD) > 1e-15,
         "An attempt to calculate a local matrix for a singular triangle (detD = " + d2s(_detD, true) + ")");

  switch (_dofs.size())
  {
  case n_dofs_first:
  {
    const double mat[][n_dofs_first] = { { 2., 1., 1. },
                                         { 1., 2., 1. },
                                         { 1., 1., 2. }
                                       };
    for (int i = 0; i < n_dofs_first; ++i)
    {
      loc_vec[i] = 0.;
      for (int j = 0; j < n_dofs_first; ++j)
        loc_vec[i] += mat[i][j] * rhs_func(points[_dofs[j]], time);
      loc_vec[i] *= fabs(_detD) / 24.;
    }
    return;
  }
  default:
    require(false, "Unknown fe order for generating local matrix");
  }
}



void Triangle::con_edge(unsigned int serial_num, unsigned int num_edge)
{
  expect(serial_num < n_edges, "");
  _edges[serial_num] = num_edge;
}



void Triangle::dis_edge(unsigned int serial_num, unsigned int num_edge)
{
  expect(serial_num < n_edges, "");
  _dis_edges[serial_num] = num_edge;
}



void Triangle::local_dg_boundary_matrix(const Edge &edge,
                                        const double coef_a,
                                        const double gamma,
                                        const DoFHandler &dof_handler,
                                        double **loc_mat) const
{
  const unsigned int dim_bound = Edge::n_vertices;
  for (int i = 0; i < dim_bound; ++i)
    for (int j = 0; j < dim_bound; ++j)
      loc_mat[i][j] = 0.0; // zero initialization

  const Point normal = edge.normal(*this, dof_handler.dofs()); // outward unit normal vector to the edge

  // we need to take a gradient of basis function associated with the first vertex.
  // so we need to find out what is its number
  Point grad_bf[Edge::n_vertices];
  for (int v = 0; v < Edge::n_vertices; ++v)
  {
    int bf = -1;
    for (int i = 0; i < n_dofs_first && bf == -1; ++i)
      if (_dofs[i] == edge.vertex(v))
        bf = i;
    expect(bf != -1, "");
    grad_bf[v] = gradient_basis_function(bf); // calculate the gradient of that basis function
    //integrals[v] = integral_over_edge(bf, edge); // calculate the integral of that basis function over the edge
  }

  // integrals of any single basis function over an edge are always
  // equal to 0.5*length_of_that_edge
  const double edge_length = edge.length(dof_handler.dofs());
  double integrals[Edge::n_vertices] = { 0.5 * edge_length,
                                         0.5 * edge_length };

  for (int i = 0; i < dim_bound; ++i)
  {
    for (int j = 0; j < dim_bound; ++j)
    {
      // B1
      loc_mat[i][j] -= coef_a * dot_product(normal, grad_bf[j]) * integrals[i];
      // B2
      loc_mat[i][j] -= coef_a * dot_product(normal, grad_bf[i]) * integrals[j];
    }
  }

  // B3
  double mat[][Edge::n_vertices] = { { 2, 1 },
                                     { 1, 2 } };
  for (int i = 0; i < dim_bound; ++i)
    for (int j = 0; j < dim_bound; ++j)
      loc_mat[i][j] += coef_a * gamma * mat[i][j] / 6.0;
}



void Triangle::local_dg_interior_matrix(const Triangle &tri,
                                        Edge edges[],
                                        double coefs[],
                                        const double gamma,
                                        const DoFHandler &dof_handler,
                                        double **loc_mat) const
{
  const unsigned int dim = Edge::n_vertices;
  for (int i = 0; i < 2 * dim; ++i)
    for (int j = 0; j < 2 * dim; ++j)
      loc_mat[i][j] = 0.0; // zero initialization

  // plus means this triangle, edge in this triangle, coefficient, etc
  // minus means another triangle 'tri', second in arrays of edges and coefs

  // auxiliary matrices
  double B_plus_plus[dim][dim];
  double B_plus_minus[dim][dim];
  double B_minus_plus[dim][dim];
  double B_minus_minus[dim][dim];

  // normals
  const Point normal_plus  = edges[0].normal(*this, dof_handler.dofs());
  const Point normal_minus = edges[1].normal(tri, dof_handler.dofs());

  expect(fabs(dot_product(normal_plus, normal_minus) + 1.0) < 1e-14,
         "Dot product of normal_plus and normal_minus is not equal to -1");

  // integrals of any single basis function over an edge are always
  // equal to 0.5*length_of_that_edge
  const double edge_length_plus  = edges[0].length(dof_handler.dofs());
  const double edge_length_minus = edges[1].length(dof_handler.dofs());
  expect(fabs(edge_length_plus - edge_length_minus) < 1e-14,
         "The lengths of dg edges corresponding to the same physical edge are not equal");
  double integrals_plus[dim]  = { 0.5 * edge_length_plus,
                                  0.5 * edge_length_plus };
  double integrals_minus[dim] = { 0.5 * edge_length_minus,
                                  0.5 * edge_length_minus };

  // we need to take gradients of basis functions associated with the vertices of the edges.
  // so we need to find out what is its number
  Point grad_bf_plus[dim];
  Point grad_bf_minus[dim];
  for (int v = 0; v < Edge::n_vertices; ++v)
  {
    // "plus"
    int bf = -1;
    for (int i = 0; i < n_dofs_first && bf == -1; ++i)
      if (this->_dofs[i] == edges[0].vertex(v))
        bf = i;
    expect(bf != -1, "");
    grad_bf_plus[v] = this->gradient_basis_function(bf); // calculate the gradient of that basis function
    // "minus"
    bf = -1;
    for (int i = 0; i < n_dofs_first && bf == -1; ++i)
      if (tri.dof(i) == edges[1].vertex(v))
        bf = i;
    expect(bf != -1, "");
    grad_bf_minus[v] = tri.gradient_basis_function(bf); // calculate the gradient of that basis function
  }

  for (int i = 0; i < dim; ++i)
  {
    for (int j = 0; j < dim; ++j)
    {
      // B1
      B_plus_plus[i][j]   = -0.5 * coefs[0] * dot_product(normal_plus,  grad_bf_plus[j])  * integrals_plus[i];
      B_plus_minus[i][j]  = -0.5 * coefs[1] * dot_product(normal_plus,  grad_bf_minus[j]) * integrals_plus[i];
      B_minus_plus[i][j]  = -0.5 * coefs[0] * dot_product(normal_minus, grad_bf_plus[j])  * integrals_minus[i];
      B_minus_minus[i][j] = -0.5 * coefs[1] * dot_product(normal_minus, grad_bf_minus[j]) * integrals_minus[i];
      // B2
      B_plus_plus[i][j]   -= 0.5 * coefs[0] * dot_product(normal_plus,  grad_bf_plus[i])  * integrals_plus[j];
      B_plus_minus[i][j]  -= 0.5 * coefs[0] * dot_product(normal_minus, grad_bf_plus[i])  * integrals_minus[j];
      B_minus_plus[i][j]  -= 0.5 * coefs[1] * dot_product(normal_plus,  grad_bf_minus[i]) * integrals_plus[j];
      B_minus_minus[i][j] -= 0.5 * coefs[1] * dot_product(normal_minus, grad_bf_minus[i]) * integrals_minus[j];
    }
  }

  // B3
  const double coef = 0.5 * (coefs[0] + coefs[1]); // average value
  double mat[][dim] = { { 2, 1 },
                        { 1, 2 } };
  for (int i = 0; i < dim; ++i)
  {
    for (int j = 0; j < dim; ++j)
    {
      B_plus_plus[i][j]   += coef * gamma * mat[i][j] / 6.0;
      B_plus_minus[i][j]  -= coef * gamma * mat[i][j] / 6.0;
      B_minus_plus[i][j]  -= coef * gamma * mat[i][j] / 6.0;
      B_minus_minus[i][j] += coef * gamma * mat[i][j] / 6.0;
    }
  }

  // assmebling local matrix
  for (int i = 0; i < dim; ++i)
  {
    for (int j = 0; j < dim; ++j)
    {
      loc_mat[i][j]         = B_plus_plus[i][j];
      loc_mat[i][j + 2]     = B_plus_minus[i][j];
      loc_mat[i + 2][j]     = B_minus_plus[i][j];
      loc_mat[i + 2][j + 2] = B_minus_minus[i][j];
    }
  }
}



const Point Triangle::gradient_basis_function(unsigned int num_bf) const
{
  expect(fabs(_detD) > 1e-15, "The triangle is singular or the barycentric coordinates are not initialized");
  // gradient of i-th basis function is (Ai, Bi, 0), since it's 2D case
  return Point(_A[num_bf], _B[num_bf], 0);
}



bool Triangle::contains_dof(unsigned int num) const
{
  for (int i = 0; i < _dofs.size(); ++i)
    if (_dofs[i] == num)
      return true;
  return false;
}
