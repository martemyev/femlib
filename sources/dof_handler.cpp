#include "dof_handler.h"
#include "fine_mesh.h"
#include "auxiliary_functions.h"
#include "sym_csr_pattern.h"
#include "finite_element.h"
#include <algorithm>


DoFHandler::DoFHandler(FineMesh *fmesh)
  : _fmesh(fmesh),
    _vertices_dofs(0)
{
  require(!_fmesh->empty(), "The fine mesh is empty");
}



DoFHandler::~DoFHandler()
{
  _cg_edges.clear();
  _discon_edges.clear();
  if (_vertices_dofs != 0)
  {
    for (unsigned i = 0; i < _fmesh->n_vertices(); ++i)
      _vertices_dofs[i].clear();
    delete[] _vertices_dofs;
  }
}



void DoFHandler::distribute_dofs(const FiniteElement &fe, COUPLING coupling)
{
  switch (fe.order())
  {
  case 1:
    {
      switch (coupling)
      {
      case CG:
        distribute_cg_first(fe);
        return;
      case DG:
        distribute_dg_first(fe);
        return;
      default:
        require(false, "Unknown coupling scheme (not CG, nor DG)");
      }
      return;
    }
  default:
    require(false, "Incorrect order of fe basis functions (" + d2s(fe.order()) + ")");
  }
}



void DoFHandler::distribute_cg_first(const FiniteElement &fe)
{
  require(fe.order() == 1, "This is not implemented for fe order other than 1");

  _dofs.resize(_fmesh->n_vertices());
  for (unsigned ver = 0; ver < _fmesh->n_vertices(); ++ver)
    _dofs[ver] = _fmesh->vertex(ver);

  // trying to work with triangles first
  for (unsigned cell = 0; cell < _fmesh->n_triangles(); ++cell)
  {
    Triangle *triangle = _fmesh->triangle_orig(cell);
    triangle->n_dofs(Triangle::n_dofs_first); // allocate the memory for degrees of freedom
    for (unsigned ver = 0; ver < Triangle::n_vertices; ++ver)
    {
      // set the numbers of degrees of freedom.
      // the number of the degree of freedom is the same as the number of the corresponding vertex
      triangle->dof(ver, triangle->vertex(ver));
    }
  }

  // then trying to work with rectangles
  for (unsigned cell = 0; cell < _fmesh->n_rectangles(); ++cell)
  {
    Rectangle *rectangle = _fmesh->rectangle_orig(cell);
    rectangle->n_dofs(Rectangle::n_dofs_first); // allocate the memory for degrees of freedom
    for (unsigned ver = 0; ver < Rectangle::n_vertices; ++ver)
    {
      // set the numbers of degrees of freedom.
      // the number of the degree of freedom is the same as the number of the corresponding vertex
      rectangle->dof(ver, rectangle->vertex(ver));
    }
  }
}



void DoFHandler::distribute_dg_first(const FiniteElement &fe)
{
  require(fe.order() == 1, "This is not implemented for fe order other than 1");

  // each triangle has n_dofs_per_triangle degrees of freedom
  const unsigned int n_dofs = Triangle::n_dofs_first * _fmesh->n_triangles();

  // allocate the memory
  _dofs.resize(n_dofs);

  // allocate memory for the map
  _vertices_dofs = new std::vector<unsigned int>[_fmesh->n_vertices()];

  // fill up the corresponding field in triangle array
  for (unsigned tr = 0; tr < _fmesh->n_triangles(); ++tr)
  {
    Triangle *triangle = _fmesh->triangle_orig(tr); // get an original triangle (not a copy)
    triangle->n_dofs(Triangle::n_dofs_first); // allocate the memory for degrees of freedom
    for (unsigned i = 0; i < Triangle::n_vertices; ++i)
    {
      const unsigned int ver_number = triangle->vertex(i); // number of vertex
      const unsigned int dof_number = Triangle::n_dofs_first * tr + i; // number of the dof
      _vertices_dofs[ver_number].push_back(dof_number); // add dof number to the corresponding list
      // set the dof
      _dofs[dof_number] = _fmesh->vertex(ver_number);
      // set the numbers of degrees of freedom.
      // the number of the degree of freedom is the same as the number of the corresponding vertex
      triangle->dof(i, dof_number);
    }
  }
}



void DoFHandler::numerate_edges()
{
  require(!_dofs.empty(), "Dofs are not distributed (initialized)!");

  SymCSRPattern dof_pattern;
  dof_pattern.make_sparse_format(*this);// make a pattern based on dofs connectivity to numerate edges which tie dofs
  _discon_edges.resize(dof_pattern.row(dof_pattern.order()));

  for (unsigned cell = 0; cell < _fmesh->n_triangles(); ++cell)
  {
    Triangle *triangle = _fmesh->triangle_orig(cell);
    unsigned int serial_n_of_edge = 0; // serial number of an edge in the structure of triangle

    for (unsigned vi = 0; vi < triangle->n_dofs(); ++vi)
    {
      const unsigned int dof_i = triangle->dof(vi);
      for (unsigned vj = 0; vj < triangle->n_dofs(); ++vj)
      {
        const unsigned int dof_j = triangle->dof(vj);
        if (dof_i > dof_j)
        {
          // the number of the edge
          const int num_edge = dof_pattern.find(dof_i, dof_j);
          expect(num_edge != -1, "");
          // add this number to the triangle structure
          triangle->dis_edge(serial_n_of_edge, num_edge);
          ++serial_n_of_edge;
          // add info to the array of DG edges
          _discon_edges[num_edge].vertex(0, dof_j); // begin vertex, since dof_j < dof_i
          _discon_edges[num_edge].vertex(1, dof_i); // end vertex, since dof_j < dof_i
          _discon_edges[num_edge].triangle(cell); // number of triangle which this edge belongs to
        }
      } // dof_j
    } // dof_i
    expect(serial_n_of_edge == Triangle::n_edges, "");
  } // cell

  SymCSRPattern ver_pattern;
  ver_pattern.make_sparse_format(*_fmesh); // make a pattern based on mesh vertices connectivity to numerate edges which tie mesh vertices - not dofs
  _cg_edges.resize(ver_pattern.row(ver_pattern.order()));

  for (unsigned cell = 0; cell < _fmesh->n_triangles(); ++cell)
  {
    Triangle *triangle = _fmesh->triangle_orig(cell);
    unsigned serial_n_of_edge = 0; // serial number of an edge in the structure of triangle

    for (unsigned vi = 0; vi < Triangle::n_vertices; ++vi)
    {
      const unsigned int ver_i = triangle->vertex(vi);
      for (unsigned vj = 0; vj < Triangle::n_vertices; ++vj)
      {
        const unsigned int ver_j = triangle->vertex(vj);
        if (ver_i > ver_j)
        {
          // the number of the edge
          const int num_edge = ver_pattern.find(ver_i, ver_j);
          expect(num_edge != -1, "");
          // add this number to the triangle structure
          triangle->con_edge(serial_n_of_edge, num_edge);
          ++serial_n_of_edge;
          // add info to the array of CG edges
          _cg_edges[num_edge].vertex(0, ver_j); // begin vertex, since ver_j < ver_i
          _cg_edges[num_edge].vertex(1, ver_i); // end vertex, since ver_j < ver_i
           // add the numbers of the DG edges that have the same location as this CG edge does
          std::vector<unsigned int> assoc_dg_edges;
          associated_dg_edges(ver_i, ver_j, _vertices_dofs, dof_pattern, assoc_dg_edges);
          _cg_edges[num_edge].edges(assoc_dg_edges);
          // change the order of some vertices in DG edges in such a way that they are oriented as CG edges are
          for (unsigned ade = 0; ade < assoc_dg_edges.size(); ++ade)
          {
            const unsigned int dg_edge = assoc_dg_edges[ade]; // the number of the associated DG edge
            // attempt to find the number of the first (begin) vertex of the DG edge
            // in the list of DG dofs of the first (begin) vertex of the CG edge
            if (find(_vertices_dofs[ver_j].begin(), _vertices_dofs[ver_j].end(), _discon_edges[dg_edge].vertex(0)) == _vertices_dofs[ver_j].end())
              _discon_edges[dg_edge].swap_vertices(); // if we couldn't find it, then we need to swap vertices of DG edge to make it be oriented as CG edge
          }
        }
      } // ver_j
    } // ver_i
    expect(serial_n_of_edge == Triangle::n_edges, "");
  } // cell
}



unsigned int DoFHandler::n_dofs() const
{
  return _dofs.size();
}



Point DoFHandler::dof(unsigned int number) const
{
  expect(number >= 0 && number < _dofs.size(), "Incorrect input parameter");
  return _dofs[number];
}



const std::vector<Point>& DoFHandler::dofs() const
{
  return _dofs;
}



const FineMesh& DoFHandler::fmesh() const
{
  return *(_fmesh);
}



unsigned int DoFHandler::vertices_dofs(unsigned int ver, unsigned int num) const
{
  expect(ver >= 0 && ver < _fmesh->n_vertices(), "");
  expect(num >= 0 && num < _vertices_dofs[ver].size(), "");
  return _vertices_dofs[ver].at(num);
}



const std::vector<unsigned int>& DoFHandler::vertices_dofs(unsigned int ver) const
{
  expect(ver >= 0 && ver < _fmesh->n_vertices(), "");
  return _vertices_dofs[ver];
}



unsigned int DoFHandler::n_dis_edges() const
{
  return _discon_edges.size();
}



unsigned int DoFHandler::n_con_edges() const
{
  return _cg_edges.size();
}



Edge DoFHandler::dis_edge(unsigned int num) const
{
  expect(num < _discon_edges.size(), "");
  return _discon_edges[num];
}



Edge DoFHandler::con_edge(unsigned int num) const
{
  expect(num < _cg_edges.size(), "");
  return _cg_edges[num];
}



void DoFHandler::boundary_dofs(const std::vector<unsigned int> &b_vertices,
                               std::vector<int> &bound_dofs) const
{
  for (unsigned bound_ver = 0; bound_ver < b_vertices.size(); ++bound_ver)
  {
    for (unsigned cell = 0; cell < _fmesh->n_triangles(); ++cell)
    {
      const Triangle tri = _fmesh->triangle(cell);
      for (unsigned ver = 0; ver < Triangle::n_vertices; ++ver)
      {
        if (tri.vertex(ver) == b_vertices[bound_ver])
          bound_dofs.push_back(tri.dof(ver));
      }
    }
  }
}



// =======================================
//
//
//
// =======================================
void associated_dg_edges(unsigned int ver_i,
                         unsigned int ver_j,
                         const std::vector<unsigned int> *vertices_dofs,
                         const SymCSRPattern &dof_pattern,
                         std::vector<unsigned int> &assoc_dg_edges)
{
  const unsigned int max_n_found_edges = 2; // for 2D case the max number of DG edges associated with CG edge equals to 2, since only 2 triangles can share an edge
  unsigned n_found_edges = 0; // current number of found edges
  for (unsigned i = 0; i < vertices_dofs[ver_i].size() && n_found_edges < max_n_found_edges; ++i)
  {
    for (unsigned j = 0; j < vertices_dofs[ver_j].size() && n_found_edges < max_n_found_edges; ++j)
    {
      const unsigned int v0 = std::min(vertices_dofs[ver_i][i], vertices_dofs[ver_j][j]); // begin vertex of the edge
      const unsigned int v1 = std::max(vertices_dofs[ver_i][i], vertices_dofs[ver_j][j]); // end vertex of the edge
      const int num_edge = dof_pattern.find(v1, v0); // attempt to find the edge
      if (num_edge != -1) // if such edge exists
      {
        assoc_dg_edges.push_back(num_edge); // add it to the list
        ++n_found_edges;
      }
    }
  }
  // the number of found edges has to be 1 (for boundary edges) or 2 (for internal edges)
  expect(n_found_edges == 1 || n_found_edges == 2, "n_found_edges is not right");
}
