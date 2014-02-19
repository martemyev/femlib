#include "mesh_element.h"
#include "auxiliary_functions.h"


MeshElement::MeshElement(unsigned int n_ver,
                         unsigned int n_edg,
                         unsigned int n_fac,
                         unsigned int el_type)
  : _n_vertices(n_ver),
    _n_edges(n_edg),
    _n_faces(n_fac),
    _gmsh_el_type(el_type),
    _material_id(0),
    _partition_id(0)
{
  _vertices.resize(_n_vertices, 0);
  _edges.resize(_n_edges, 0);
  _faces.resize(_n_faces, 0);
}



MeshElement::~MeshElement()
{
  _vertices.clear();
  _edges.clear();
  _faces.clear();
}



unsigned int MeshElement::n_vertices() const
{
  expect(_n_vertices == _vertices.size(),
         "Memory for vertices is not allocated properly (size is " + d2s(_vertices.size()) +
         "), or n_vertices (" + d2s(_n_vertices) + ") is set to wrong number");
  return _n_vertices;
}



unsigned int MeshElement::n_edges() const
{
  expect(_n_edges == _edges.size(),
         "Memory for edges is not allocated properly (size is " + d2s(_edges.size()) +
         "), or n_edges (" + d2s(_n_edges) + ") is set to wrong number");
  return _n_edges;
}



unsigned int MeshElement::n_faces() const
{
  expect(_n_faces == _faces.size(),
         "Memory for faces is not allocated properly (size is " + d2s(_faces.size()) +
         "), or n_faces (" + d2s(_n_faces) + ") is set to wrong number");
  return _n_faces;
}



unsigned int MeshElement::gmsh_el_type() const
{
  return _gmsh_el_type;
}



unsigned int MeshElement::material_id() const
{
  return _material_id;
}



unsigned int MeshElement::partition_id() const
{
  return _partition_id;
}



MeshElement::MeshElement(const MeshElement &elem)
  : _n_vertices(elem._n_vertices),
    _n_edges(elem._n_edges),
    _n_faces(elem._n_faces),
    _gmsh_el_type(elem._gmsh_el_type),
    _material_id(elem._material_id),
    _partition_id(elem._partition_id)
{
  _vertices = elem._vertices;
  _edges = elem._edges;
  _faces = elem._faces;
}



MeshElement& MeshElement::operator =(const MeshElement &elem)
{
  _n_vertices = elem._n_vertices;
  _n_edges = elem._n_edges;
  _n_faces = elem._n_faces;
  _gmsh_el_type = elem._gmsh_el_type;
  _material_id = elem._material_id;
  _partition_id = elem._partition_id;
  _vertices = elem._vertices;
  _edges = elem._edges;
  _faces = elem._faces;
  return *this;
}



unsigned int MeshElement::vertex(unsigned int number) const
{
  expect(number < _n_vertices,
         "The local number of vertex is incorrect: " + d2s(number) +
         ". It has to be in range [0, " + d2s(_n_vertices) + ").");
  return _vertices[number];
}



unsigned int MeshElement::edge(unsigned int number) const
{
  expect(number < _n_edges,
         "The local number of edge is incorrect: " + d2s(number) +
         ". It has to be in range [0, " + d2s(_n_edges) + ").");
  return _edges[number];
}



unsigned int MeshElement::face(unsigned int number) const
{
  expect(number < _n_faces,
         "The local number of face is incorrect: " + d2s(number) +
         ". It has to be in range [0, " + d2s(_n_faces) + ").");
  return _faces[number];
}



void MeshElement::vertex(unsigned int local_number, unsigned int global_number)
{
  expect(local_number < n_vertices(),
         "Local number (" + d2s(local_number) +
         ") is incorrect. It must be in the range [0, " + d2s(_n_vertices) + ")");
  _vertices[local_number] = global_number;
}



void MeshElement::edge(unsigned int local_number, unsigned int global_number)
{
  expect(local_number < n_edges(),
         "Local number (" + d2s(local_number) +
         ") is incorrect. It must be in the range [0, " + d2s(_n_edges) + ")");
  _edges[local_number] = global_number;
}



void MeshElement::face(unsigned int local_number, unsigned int global_number)
{
  expect(local_number < n_faces(),
         "Local number (" + d2s(local_number) +
         ") is incorrect. It must be in the range [0, " + d2s(_n_faces) + ")");
  _faces[local_number] = global_number;
}



void MeshElement::faces(const std::vector<unsigned int> &face_numbers)
{
  expect(face_numbers.size() == n_faces(),
         "Array of face numbers has another size (" + d2s(face_numbers.size()) +
         ") than it must be (" + d2s(n_faces()) + ")");
  _faces = face_numbers;
}



bool MeshElement::contains(const unsigned int vertex) const
{
  for (unsigned int i = 0; i < _n_vertices; ++i)
    if (vertex == _vertices[i])
      return true;
  return false;
}

