#include "line.h"
#include "auxiliary_functions.h"

NAMESPACE_FEM_OPEN



Line::Line()
  : MeshElement(n_vertices, n_edges, n_faces, gmsh_el_type)
{ }



Line::Line(const std::vector<unsigned int> &ver,
           const unsigned int mat_id)
  : MeshElement(n_vertices, n_edges, n_faces, gmsh_el_type)
{
  expect(_vertices.size() == ver.size(),
         "The size of list of vertices (" + d2s(ver.size()) +
         ") is not equal to really needed number of vertices (" + d2s(n_vertices) + ")");
  _vertices = ver;
  _material_id = mat_id;
}



Line::Line(const unsigned int v1,
           const unsigned int v2,
           const unsigned int mat_id)
  : MeshElement(n_vertices, n_edges, n_faces, gmsh_el_type)
{
  _vertices[0] = v1;
  _vertices[1] = v2;
  _material_id = mat_id;
}



Line::~Line()
{ }



unsigned int Line::common_vertex(const Line& line) const
{
  for (unsigned int i = 0; i < n_vertices; ++i)
    if (line.contains(_vertices[i]))
      return _vertices[i];
  require(false, "There is no common vertex between these two lines!");
  return 0; // to calm compiler down
}



unsigned int Line::another_vertex(const unsigned int vertex) const
{
  if (vertex == _vertices[0])
    return _vertices[1];
  else if (vertex == _vertices[1])
    return _vertices[0];
  else
    require(false, "This line doesn't contain the vertex. So we can't find another one.");
  return 0; // to calm compiler down
}


NAMESPACE_FEM_CLOSE
