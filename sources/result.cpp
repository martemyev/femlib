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

  const std::vector<Point> dofs = _dof_handler->dofs(); // the list of all degrees of freedom
  const FineMesh *fmesh = _dof_handler->fmesh(); // the fine triangular mesh

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
  out << "    <Piece NumberOfPoints=\"" << dofs.size() << "\" NumberOfCells=\"" << fmesh->n_triangles() << "\">\n";
  out << "      <PointData>\n";
  out << "        <DataArray type=\"Float64\" Name=\"solution\" format=\"ascii\">\n";
  for (int i = 0; i < dofs.size(); ++i)
    out << solution_values[i] << "\n";
  out << "        </DataArray>\n";
  if (exact_solution)
  {
    out << "        <DataArray type=\"Float64\" Name=\"exact_solution\" format=\"ascii\">\n";
    for (int i = 0; i < dofs.size(); ++i)
      out << exact_solution_values[i] << "\n";
    out << "        </DataArray>\n";
  }
  out << "      </PointData>\n";
  out << "      <Points>\n";
  out << "        <DataArray type=\"Float64\" NumberOfComponents=\"" << Point::n_coord << "\" format=\"ascii\">\n";
  for (int i = 0; i < dofs.size(); ++i)
  {
    for (int j = 0; j < Point::n_coord; ++j)
      out << dofs[i].coord(j) << " ";
    out << "\n";
  }
  out << "        </DataArray>\n";
  out << "      </Points>\n";
  out << "      <Cells>\n";
  out << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
  for (int tr = 0; tr < fmesh->n_triangles(); ++tr)
  {
    const Triangle triangle = fmesh->triangle(tr);
    for (int d = 0; d < Triangle::n_dofs_first; ++d)
      out << triangle.dof(d) << " ";
    out << "\n";
  }
  out << "        </DataArray>\n";
  out << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  for (int tr = 0; tr < fmesh->n_triangles(); ++tr)
    out << (tr + 1) * Triangle::n_dofs_first << " ";
  out << "\n";
  out << "        </DataArray>\n";
  out << "        <DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">\n";
  for (int tr = 0; tr < fmesh->n_triangles(); ++tr)
    out << Triangle::vtk_el_type << " ";
  out << "\n";
  out << "        </DataArray>\n";
  out << "      </Cells>\n";
  out << "    </Piece>\n";
  out << "  </UnstructuredGrid>\n";
  out << "</VTKFile>\n";

  out.close();
}