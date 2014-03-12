#include "result.h"
#include "auxiliary_functions.h"
#include "dof_handler.h"
#include "fine_mesh.h"
#include "point.h"
#include <boost/filesystem.hpp>
#include <fstream>



Result::Result(const DoFHandler *dof_handler)
  : _dof_handler(dof_handler)
{ }



Result::~Result()
{ }



void Result::write_vtu(const std::string &filename,
                       const Vec &solution,
                       const Vec &exact_solution) const
{
  using namespace boost::filesystem;
  expect(extension(filename) == ".vtu",
         "The extension of the file ('" + extension(filename) +
         "') is not suitable for this function, because it should be '.vtu'.");

//  expect(_dof_handler->n_dofs() == _dof_handler->fmesh()->n_vertices(),
//         "This function should be rewritten for the case higher order basis functions"
//         "(when the number of degrees of freedom is not equal to the number of the mesh vertices)");

  std::ofstream out(filename.c_str()); // open the file for writing
  require(out, "File " + filename + " cannot be opened");

  const std::vector<Point> &dofs = _dof_handler->dofs(); // the list of all degrees of freedom
  require(!dofs.empty(), "The vector of degrees of freedom is empty");

  const FineMesh &fmesh = _dof_handler->fmesh(); // the fine triangular mesh
  //require(fmesh != 0, "Fine mesh is not initialized for this dof handler");

  // extract the data from PETSc vectors
  std::vector<int> idx(dofs.size());
  std::iota(idx.begin(), idx.end(), 0); // idx = { 0, 1, 2, 3, .... }
  std::vector<double> solution_values(dofs.size());
  VecGetValues(solution, dofs.size(), &idx[0], &solution_values[0]); // extract the values of the numerical solution

  std::vector<double> exact_solution_values;
  if (exact_solution)
  {
    exact_solution_values.resize(dofs.size());
    VecGetValues(exact_solution, dofs.size(), &idx[0], &exact_solution_values[0]); // extract the values of the exact solution
  }

  out << "<?xml version=\"1.0\"?>\n";
  out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  out << "  <UnstructuredGrid>\n";
  out << "    <Piece NumberOfPoints=\"" << dofs.size() << "\" NumberOfCells=\"" << fmesh.n_triangles() << "\">\n";
  out << "      <PointData>\n";
  out << "        <DataArray type=\"Float64\" Name=\"solution\" format=\"ascii\">\n";
  for (unsigned int i = 0; i < dofs.size(); ++i)
    out << solution_values[i] << "\n";
  out << "        </DataArray>\n";
  if (exact_solution)
  {
    out << "        <DataArray type=\"Float64\" Name=\"exact_solution\" format=\"ascii\">\n";
    for (unsigned int i = 0; i < dofs.size(); ++i)
      out << exact_solution_values[i] << "\n";
    out << "        </DataArray>\n";
  }
  out << "      </PointData>\n";
  out << "      <Points>\n";
  out << "        <DataArray type=\"Float64\" NumberOfComponents=\"" << Point::n_coord << "\" format=\"ascii\">\n";
  for (unsigned int i = 0; i < dofs.size(); ++i)
  {
    for (unsigned int j = 0; j < Point::n_coord; ++j)
      out << dofs[i].coord(j) << " ";
    out << "\n";
  }
  out << "        </DataArray>\n";
  out << "      </Points>\n";
  out << "      <Cells>\n";
  out << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
  for (unsigned int tr = 0; tr < fmesh.n_triangles(); ++tr)
  {
    const Triangle triangle = fmesh.triangle(tr);
    for (unsigned int d = 0; d < Triangle::n_dofs_first; ++d)
      out << triangle.dof(d) << " ";
    out << "\n";
  }
  out << "        </DataArray>\n";
  out << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  for (unsigned int tr = 0; tr < fmesh.n_triangles(); ++tr)
    out << (tr + 1) * Triangle::n_dofs_first << " ";
  out << "\n";
  out << "        </DataArray>\n";
  out << "        <DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">\n";
  for (unsigned int tr = 0; tr < fmesh.n_triangles(); ++tr)
    out << Triangle::vtk_el_type << " ";
  out << "\n";
  out << "        </DataArray>\n";
  out << "      </Cells>\n";
  out << "    </Piece>\n";
  out << "  </UnstructuredGrid>\n";
  out << "</VTKFile>\n";

  out.close();
}



void Result::write_vts(const std::string &filename,
                       unsigned int N_FINE_X,
                       unsigned int N_FINE_Y,
                       const Vec &solution,
                       const Vec &exact_solution,
                       const std::vector<double> &coef_alpha,
                       const std::vector<double> &coef_beta) const
{
  using namespace boost::filesystem;
  expect(extension(filename) == ".vts",
         "The extension of the file ('" + extension(filename) +
         "') is not suitable for this function, because it should be '.vts'.");

//  expect(_dof_handler->n_dofs() == _dof_handler->fmesh()->n_vertices(),
//         "This function should be rewritten for the case higher order basis functions"
//         "(when the number of degrees of freedom is not equal to the number of the mesh vertices)");

  std::ofstream out(filename.c_str()); // open the file for writing
  require(out, "File " + filename + " cannot be opened");

  const std::vector<Point> &dofs = _dof_handler->dofs(); // the list of all degrees of freedom
  require(!dofs.empty(), "The vector of degrees of freedom is empty");

  const FineMesh &fmesh = _dof_handler->fmesh(); // the fine triangular mesh
  //require(fmesh != 0, "Fine mesh is not initialized for this dof handler");

  // extract the data from PETSc vectors
  std::vector<int> idx(dofs.size());
  std::iota(idx.begin(), idx.end(), 0); // idx = { 0, 1, 2, 3, .... }
  std::vector<double> solution_values(dofs.size());
  VecGetValues(solution, dofs.size(), &idx[0], &solution_values[0]); // extract the values of the numerical solution

  std::vector<double> exact_solution_values;
  if (exact_solution)
  {
    exact_solution_values.resize(dofs.size());
    VecGetValues(exact_solution, dofs.size(), &idx[0], &exact_solution_values[0]); // extract the values of the exact solution
  }

  out << "<?xml version=\"1.0\"?>\n";
  out << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  out << "  <StructuredGrid WholeExtent=\"1 " << N_FINE_X + 1 << " 1 " << N_FINE_Y + 1 << " 1 1\">\n";
  out << "    <Piece Extent=\"1 " << N_FINE_X + 1 << " 1 " << N_FINE_Y + 1 << " 1 1\">\n";
  out << "      <PointData Scalars=\"scalars\">\n";
  out << "        <DataArray type=\"Float64\" Name=\"U_solution\" format=\"ascii\">\n";
  for (unsigned int i = 0; i < dofs.size(); ++i)
    out << solution_values[i] << "\n";
  out << "        </DataArray>\n";
  if (exact_solution)
  {
    out << "        <DataArray type=\"Float64\" Name=\"U_exact\" format=\"ascii\">\n";
    for (unsigned int i = 0; i < dofs.size(); ++i)
      out << exact_solution_values[i] << "\n";
    out << "        </DataArray>\n";
  }
  out << "      </PointData>\n";
  if (!coef_alpha.empty() || !coef_beta.empty())
  {
    out << "      <CellData Scalars=\"scalars\">\n";
    if (!coef_alpha.empty())
    {
      out << "        <DataArray type=\"Float64\" Name=\"coef_alpha\" format=\"ascii\">\n";
      expect(coef_alpha.size() == N_FINE_X * N_FINE_Y, "dimensions mismatch. look at the code");
      for (unsigned int i = 0; i < coef_alpha.size(); ++i)
        out << coef_alpha[i] << "\n";
      out << "        </DataArray>\n";
    }
    if (!coef_beta.empty())
    {
      out << "        <DataArray type=\"Float64\" Name=\"coef_beta\" format=\"ascii\">\n";
      expect(coef_beta.size() == N_FINE_X * N_FINE_Y, "dimensions mismatch. look at the code");
      for (unsigned int i = 0; i < coef_beta.size(); ++i)
        out << coef_beta[i] << "\n";
      out << "        </DataArray>\n";
    }
    out << "      </CellData>\n";
  }
  out << "      <Points>\n";
  out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (unsigned int i = 0; i < fmesh.n_vertices(); ++i)
  {
    const Point vert = fmesh.vertex(i);
    out << vert.coord(0) << " " << vert.coord(1) << " " << vert.coord(2) << " ";
  }
  out << "\n";
  out << "        </DataArray>\n";
  out << "      </Points>\n";
  out << "    </Piece>\n";
  out << "  </StructuredGrid>\n";
  out << "</VTKFile>\n";

  out.close();
}

