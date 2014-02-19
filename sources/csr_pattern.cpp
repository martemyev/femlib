#include "csr_pattern.h"
#include "triangle.h"
#include "fine_mesh.h"
#include "auxiliary_functions.h"


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

  // the number of rows of the matrix (and its order) of connectivity between degrees of freedom
  _order = dof_handler.n_dofs();

  std::set<unsigned int> *connect = new std::set<unsigned int>[_order];

  // pass through all triangles and all dofs on them
  expect(dof_handler.fmesh()->n_triangles() != 0, "");
  for (int cell = 0; cell < dof_handler.fmesh()->n_triangles(); ++cell)
  {
    Triangle triangle = dof_handler.fmesh()->triangle(cell);

    expect(triangle.n_dofs() != 0, "");
    for (int di = 0; di < triangle.n_dofs(); ++di)
    {
      const unsigned int dof_i = triangle.dof(di); // the number of the first degree of freedom
      for (int dj = 0; dj < triangle.n_dofs(); ++dj)
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
  for (int i = 0; i < _order; ++i)
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
  expect(dof_handler.fmesh()->n_triangles() != 0, "");
  for (int cell = 0; cell < dof_handler.fmesh()->n_triangles(); ++cell)
  {
    Triangle triangle = dof_handler.fmesh()->triangle(cell);

    expect(triangle.n_dofs() != 0, "");
    for (int di = 0; di < triangle.n_dofs(); ++di)
    {
      const unsigned int dof_i = triangle.dof(di); // the number of the first degree of freedom
      for (int dj = 0; dj < triangle.n_dofs(); ++dj)
      {
        const unsigned int dof_j = triangle.dof(dj); // the number of the second degree of freedom
        // insert the values in the corresponding places
        connect[dof_i].insert(dof_j);
        connect[dof_j].insert(dof_i);
      }
    }
  }

  // pass through all interior CG edges
  for (int e = 0; e < dof_handler.n_con_edges(); ++e)
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
  for (int i = 0; i < _order; ++i)
    connect[i].clear();
  delete[] connect;
}




void CSRPattern::pattern_initialization(const std::set<unsigned int> *connect)
{
  require(_row.empty(), "Row vector is not empty");
  require(_col.empty(), "Col vector is not empty");

  _row.resize(_order + 1);
  _row[0] = 0;
  for (int i = 0; i < _order; ++i)
    _row[i + 1] = _row[i] + connect[i].size();

  _col.resize(_row[_order]);
  int k = 0;
  for (int i = 0; i < _order; ++i)
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
  for (int i = _row[num_row]; i < _row[num_row + 1]; ++i)
    if (num_col == _col[i])
      return i;

  //require(false, "The serial number of nonzero element cannot be found");
  //return 0; // to calm down a compiler
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
  for (int i = 0; i < _order; ++i)
    nnz[i] = _row[i + 1] - _row[i];
  return nnz;
}
