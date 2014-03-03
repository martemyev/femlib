#include "quadrangle.h"


Quadrangle::Quadrangle()
  : MeshElement(n_vertices, n_edges, n_faces, gmsh_el_type)
{ }



Quadrangle::Quadrangle(const std::vector<unsigned int> &ver,
                       const unsigned int mat_id)
  : MeshElement(n_vertices, n_edges, n_faces, gmsh_el_type)
{
  expect(vertices.size() == ver.size(),
         "The size of list of vertices (" + d2s(ver.size()) +
         ") is not equal to really needed number of vertices (" + d2s(n_vertices) + ")");
  vertices = ver;
  material_id = mat_id;
}



Quadrangle::Quadrangle(const unsigned int v1,
                       const unsigned int v2,
                       const unsigned int v3,
                       const unsigned int v4,
                       const unsigned int mat_id)
  : MeshElement(n_vertices, n_edges, n_faces, gmsh_el_type)
{
  vertices[0] = v1;
  vertices[1] = v2;
  vertices[2] = v3;
  vertices[3] = v4;
  material_id = mat_id;
}



Quadrangle::~Quadrangle()
{ }
