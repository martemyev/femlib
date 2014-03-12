#include "quadrangle.h"
#include "auxiliary_functions.h"



Quadrangle::Quadrangle()
  : MeshElement(n_vertices, n_edges, n_faces, gmsh_el_type)
{
  for (unsigned int i = 0; i < n_vertices; ++i)
    _X[i] = _Y[i] = 0.;
}



Quadrangle::~Quadrangle()
{
  //_ghost_cells.clear();
  //_dofs.clear();
  //_dis_edges.clear();
}



Quadrangle::Quadrangle(const std::vector<unsigned int> &ver,
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
  //_ghost_cells = g_cells;
  //_dis_edges.resize(n_edges, 0);

  if (!mesh_vertices.empty())
  {
    for (unsigned int i = 0; i < n_vertices; ++i)
    {
      expect(_vertices[i] < mesh_vertices.size(), ""); // we can check this way, since the vertices in mesh_vertices vector have dense numeration from 0
      _X[i] = mesh_vertices[_vertices[i]].coord(0);
      _Y[i] = mesh_vertices[_vertices[i]].coord(1);
    }
  }
}



Quadrangle::Quadrangle(const unsigned int v1,
                       const unsigned int v2,
                       const unsigned int v3,
                       const unsigned int v4,
                       const std::vector<Point> &mesh_vertices,
                       const unsigned int mat_id,
                       const unsigned int part_id,
                       const std::vector<unsigned int> &g_cells)
  : MeshElement(n_vertices, n_edges, n_faces, gmsh_el_type)
{
  _vertices[0] = v1;
  _vertices[1] = v2;
  _vertices[2] = v3;
  _vertices[3] = v4;
  _material_id = mat_id;
  _partition_id = part_id;
  //_ghost_cells = g_cells;
  //_dis_edges.resize(n_edges, 0);

  if (!mesh_vertices.empty())
  {
    for (unsigned int i = 0; i < n_vertices; ++i)
    {
      expect(_vertices[i] < mesh_vertices.size(), ""); // we can check this way, since the vertices in mesh_vertices vector have dense numeration from 0
      _X[i] = mesh_vertices[_vertices[i]].coord(0);
      _Y[i] = mesh_vertices[_vertices[i]].coord(1);
    }
  }
}



Quadrangle::Quadrangle(const Quadrangle &quad)
  : MeshElement(quad)//,
    //_ghost_cells(tri._ghost_cells),
    //_dofs(tri._dofs),
    //_detD(tri._detD)
{
  for (unsigned int i = 0; i < n_vertices; ++i)
  {
    _X[i] = quad._X[i];
    _Y[i] = quad._Y[i];
  }
  //_dis_edges = tri._dis_edges;
}



Quadrangle& Quadrangle::operator =(const Quadrangle &quad)
{
  MeshElement::operator =(quad);
  //_ghost_cells = tri._ghost_cells;
  //_dofs = tri._dofs;
  //_detD = tri._detD;
  for (unsigned int i = 0; i < n_vertices; ++i)
  {
    _X[i] = quad._X[i];
    _Y[i] = quad._Y[i];
  }
  //_dis_edges = tri._dis_edges;

  return *this;
}
