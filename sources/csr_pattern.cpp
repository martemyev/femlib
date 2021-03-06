#include "csr_pattern.h"
#include "triangle.h"
#include "fine_mesh.h"
#include "auxiliary_functions.h"

NAMESPACE_FEM_OPEN



CSRPattern::CSRPattern()
{ }



CSRPattern::~CSRPattern()
{
  clear();
}



void CSRPattern::clear()
{
  _row.clear();
  _col.clear();
}



void CSRPattern::make_sparse_format(const DoFHandler &dof_handler, COUPLING coupling)
{
  switch (coupling)
  {
  case CG:
    make_cg_sparse_format(dof_handler);
    return;
  case DG:
    make_dg_sparse_format(dof_handler);
    return;
  default:
    require(false, "Unknown coupling (nor CG, nor DG)");
  }
}



void CSRPattern::make_cg_sparse_format(const DoFHandler &dof_handler)
{
  clear(); // in the case that sparse format was already initialized somehow

  require(dof_handler.n_dofs() != 0, "Dofs are not initialized yet");

  if (dof_handler.fmesh().n_triangles() == 0)
  {
    make_sparse_format(dof_handler.fmesh().rectangles(), dof_handler.n_dofs(), DOFS);
    return;
  }

  // the number of rows of the matrix (and its order) of connectivity between degrees of freedom
  _order = dof_handler.n_dofs();

  std::set<unsigned int> *connect = new std::set<unsigned int>[_order];

  // pass through all triangles and all dofs on them
  expect(dof_handler.fmesh().n_triangles() != 0, "");
  for (unsigned int cell = 0; cell < dof_handler.fmesh().n_triangles(); ++cell)
  {
    Triangle triangle = dof_handler.fmesh().triangle(cell);

    expect(triangle.n_dofs() != 0, "");
    for (unsigned int di = 0; di < triangle.n_dofs(); ++di)
    {
      const unsigned int dof_i = triangle.dof(di); // the number of the first degree of freedom
      for (unsigned int dj = 0; dj < triangle.n_dofs(); ++dj)
      {
        const unsigned int dof_j = triangle.dof(dj); // the number of the second degree of freedom
        // insert the values in the corresponding places
        connect[dof_i].insert(dof_j);
        connect[dof_j].insert(dof_i);
      }
    }
  }

  // initialization of the pattern
  pattern_initialization(connect);

  // free the memory
  for (unsigned int i = 0; i < _order; ++i)
    connect[i].clear();
  delete[] connect;
}



void CSRPattern::make_dg_sparse_format(const DoFHandler &dof_handler)
{
  clear(); // in the case that sparse format was already initialized somehow

  require(dof_handler.n_dofs() != 0, "Dofs are not initialized yet");
  require(dof_handler.n_con_edges() != 0, "Edges are not numerated yet");

  // the number of rows of the matrix (and its order) of connectivity between degrees of freedom
  _order = dof_handler.n_dofs();

  std::set<unsigned int> *connect = new std::set<unsigned int>[_order];

  // pass through all triangles and all dofs on them
  expect(dof_handler.fmesh().n_triangles() != 0, "");
  for (unsigned int cell = 0; cell < dof_handler.fmesh().n_triangles(); ++cell)
  {
    Triangle triangle = dof_handler.fmesh().triangle(cell);

    expect(triangle.n_dofs() != 0, "");
    for (unsigned int di = 0; di < triangle.n_dofs(); ++di)
    {
      const unsigned int dof_i = triangle.dof(di); // the number of the first degree of freedom
      for (unsigned int dj = 0; dj < triangle.n_dofs(); ++dj)
      {
        const unsigned int dof_j = triangle.dof(dj); // the number of the second degree of freedom
        // insert the values in the corresponding places
        connect[dof_i].insert(dof_j);
        connect[dof_j].insert(dof_i);
      }
    }
  }

  // pass through all interior CG edges
  for (unsigned int e = 0; e < dof_handler.n_con_edges(); ++e)
  {
    const Edge edge = dof_handler.con_edge(e); // CG edge
    if (edge.n_assoc_edges() == 2) // for interior edges only
    {
      Edge dg_edges[2];
      dg_edges[0] = dof_handler.dis_edge(edge.assoc_edge(0));
      dg_edges[1] = dof_handler.dis_edge(edge.assoc_edge(1));
      unsigned int dofs[4];
      dofs[0] = dg_edges[0].vertex(0);
      dofs[1] = dg_edges[0].vertex(1);
      dofs[2] = dg_edges[1].vertex(0);
      dofs[3] = dg_edges[1].vertex(1);
      for (int i = 0; i < 4; ++i)
      {
        for (int j = 0; j < 4; ++j)
        {
          connect[dofs[i]].insert(dofs[j]);
          connect[dofs[j]].insert(dofs[i]);
        }
      }
    } // interior edges
  } // all edges

  // initialization of the pattern
  pattern_initialization(connect);

  // free the memory
  for (unsigned int i = 0; i < _order; ++i)
    connect[i].clear();
  delete[] connect;
}



void CSRPattern::make_sparse_format(const std::vector<Rectangle> &rectangles, unsigned int order, CONNECTIVITY connectivity)
{
  clear(); // free some memory in the case if sparse format was already made

  // the number of rows of the matrix (and its order) of connectivity between vertices or degrees of freedom
  _order = order;

  std::set<unsigned int> *connect = new std::set<unsigned int>[_order];

  // pass through all rectangles and all dofs on them
  for (unsigned int cell = 0; cell < rectangles.size(); ++cell)
  {
    const Rectangle rect = rectangles[cell];

    if (connectivity == DOFS)
      expect(rect.n_dofs() != 0, "");

    unsigned int N = -1; // the local number of vertices or degrees of freedom
    if (connectivity == VERTICES)
      N = Rectangle::n_vertices;
    else if (connectivity == DOFS)
      N = rect.n_dofs();
    else
      require(false, "");

    for (unsigned int ii = 0; ii < N; ++ii)
    {
      unsigned int num_i = -1;
      if (connectivity == VERTICES)
        num_i = rect.vertex(ii); // the number of the first vertex
      else if (connectivity == DOFS)
        num_i = rect.dof(ii); // the number of the first degree of freedom
      else
        require(false, "");

      for (unsigned int jj = 0; jj < N; ++jj)
      {
        unsigned int num_j = -1;
        if (connectivity == VERTICES)
          num_j = rect.vertex(jj); // the number of the second vertex
        else if (connectivity == DOFS)
          num_j = rect.dof(jj); // the number of the second degree of freedom
        else
          require(false, "");

        // insert the values in the corresponding places
        connect[num_i].insert(num_j);
        connect[num_j].insert(num_i);
      }
    }
  }

  // initialization of the pattern
  pattern_initialization(connect);

  // free the memory
  for (unsigned int i = 0; i < _order; ++i)
    connect[i].clear();
  delete[] connect;
}




void CSRPattern::pattern_initialization(const std::set<unsigned int> *connect)
{
  require(_row.empty(), "Row vector is not empty");
  require(_col.empty(), "Col vector is not empty");

  _row.resize(_order + 1);
  _row[0] = 0;
  for (unsigned int i = 0; i < _order; ++i)
    _row[i + 1] = _row[i] + connect[i].size();

  _col.resize(_row[_order]);
  int k = 0;
  for (unsigned int i = 0; i < _order; ++i)
  {
    for (std::set<unsigned int>::const_iterator iter = connect[i].begin();
         iter != connect[i].end();
         ++iter)
    {
      _col[k] = *iter;
      ++k;
    }
  }
}



int CSRPattern::find(unsigned int num_row, unsigned num_col) const
{
  for (unsigned int i = _row[num_row]; i < _row[num_row + 1]; ++i)
    if (num_col == _col[i])
      return i;

  // if we found nothing
  return -1;
}



unsigned int CSRPattern::order() const
{
  return _order;
}



unsigned int CSRPattern::row(unsigned int number) const
{
  expect(number >= 0 && number < _row.size(), "Incorrect input data");
  return _row[number];
}



unsigned int CSRPattern::col(unsigned int number) const
{
  expect(number >= 0 && number < _col.size(), "Incorrect input data");
  return _col[number];
}



const int* CSRPattern::nnz() const
{
  int *nnz = new int[_order];
  for (unsigned int i = 0; i < _order; ++i)
    nnz[i] = _row[i + 1] - _row[i];
  return nnz;
}


NAMESPACE_FEM_CLOSE
