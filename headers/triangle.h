#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "mesh_element.h"
#include "point.h"
#include "parameters.h"

class Edge;
class DoFHandler;


/**
 * Triangle - a 2-dimensional simplex.
 * The simplest shape in the 2D space.
 * It's an element of the mesh,
 * therefore it inherits the most part of
 * functionality from the MeshElement class.
 */
class Triangle : public MeshElement
{
public:
            /**
             * The number of vertices of triangle
             */
  static const unsigned int n_vertices = 3;

            /**
             * The number of edges of triangle
             */
  static const unsigned int n_edges = 3;

            /**
             * Triangle is a 2D shape,
             * so it's a face itself (for a tetrahedron)
             */
  static const unsigned int n_faces = 1;

            /**
             * In Gmsh a triangle is defined by number 2
             */
  static const unsigned int gmsh_el_type = 2;

            /**
             * In VTK format a triangle is defined by number 5
             */
  static const unsigned int vtk_el_type = 5;

            /**
             * The number of degrees of freedom on the triangle
             * in case of first order basis functions coincides
             * with the number of vertices
             */
  static const unsigned int n_dofs_first = n_vertices;

            /**
             * Default constructor
             */
  Triangle();

            /**
             * Destructor
             */
  virtual ~Triangle();

            /**
             * Constructor with parameters
             * @param ver - triangle vertices
             * @param mat_id - material ID
             * @param part_id - partition ID (can be used as a number of a coarse element)
             * @param ghost_cells - the numbers of other partitions (coarse element)
             *                      that are connected with this triangle at least by one vertex
             */
  Triangle(const std::vector<unsigned int> &ver,
           const std::vector<Point> &mesh_vertices = std::vector<Point>(),
           const unsigned int mat_id = 0,
           const unsigned int part_id = 0,
           const std::vector<unsigned int> &ghost_cells = std::vector<unsigned int>());

            /**
             * Constructor with parameters
             * @param v1 - first vertex
             * @param v2 - second vertex
             * @param v3 - third vertex
             * @param mat_id - material ID
             * @param part_id - partition ID (can be used as a number of a coarse element)
             * @param ghost_cells - the numbers of other partitions (coarse element)
             *                      that are connected with this triangle at least by one vertex
             */
  Triangle(const unsigned int v1,
           const unsigned int v2,
           const unsigned int v3,
           const std::vector<Point> &mesh_vertices = std::vector<Point>(),
           const unsigned int mat_id = 0,
           const unsigned int part_id = 0,
           const std::vector<unsigned int> &ghost_cells = std::vector<unsigned int>());

            /**
             * Copy constructor
             */
  Triangle(const Triangle &tri);

            /**
             * Copy assignment operator
             */
  Triangle& operator =(const Triangle &tri);

            /**
             * Get the number of the ghost cells
             */
  unsigned int n_ghost_cells() const;

            /**
             * Get the number of 'num'-th ghost cell
             * @param num - the serial number of ghost cell from the vector
             */
  unsigned int ghost_cell(unsigned int num) const;

            /**
             * Get the number of degrees of freedom associated with this triangle
             */
  unsigned int n_dofs() const;

            /**
             * Get the number of degree of freedom
             * @param number - the serial number of the dof for this triangle
             */
  unsigned int dof(unsigned int number) const;

            /**
             * Set the number of degrees of freedom on the triangle.
             * Actually - just to allocate the memory for _dofs vector
             * @param number - the number of degrees of freedom associated with the triangle
             */
  void n_dofs(unsigned int number);

            /**
             * Set the value of a degree of freedom
             * @param number - the serial number of the dof
             * @param value - the value of the dof
             */
  void dof(unsigned int number, unsigned int value);

            /**
             * Generate the local mass matrix for the triangle
             * @param loc_mat - output data, generated local matrix
             */
  void local_mass_matrix(double **loc_mat, double coef_alpha) const;

            /**
             * Generate the local stiffness matrix for the triangle
             * @param loc_mat - output data, generated local matrix
             * @param coef_a - the value of the coefficient a on this triangle
             */
  void local_stiffness_matrix(double **loc_mat, double coef_beta) const;

            /**
             * Generate the local matrix of integrals over a one boundary edge belogning to this triangle
             * @param loc_mat - output data, generated local matrix
             * @param edge - the edge over which we need to take the integrals
             * @param coef_a - the coefficient a on this triangle
             * @param gamma - the IPDG parameter
             * @param dof_handler - the constant reference to an object of DoFHandler class
             *                      that provides us with the list of all degrees of freedom, for example
             */
  void local_dg_boundary_matrix(double **loc_mat,
                                const Edge &edge,
                                double coef_a,
                                double gamma,
                                const DoFHandler &dof_handler) const;

            /**
             * Generate the local matrix of integrals over a one interior CG edge sharing by this triangle and another one
             * @param loc_mat - output data, generated local matrix
             * @param tri - another triangle
             * @param edges - the DG edges: edges[0] belongs to this triangle, edges[1] - to another one
             * @param coefs - the coefficients a: coefs[0] is on this triangle, coefs[1] is on another one
             * @param gamma - the IPDG parameter
             * @param dof_handler - the constant reference to an object of DoFHandler class
             *                      that provides us with the list of all degrees of freedom, for example
             */
  void local_dg_interior_matrix(double **loc_mat,
                                const Triangle &tri,
                                Edge edges[],
                                double coefs[],
                                double gamma,
                                const DoFHandler &dof_handler) const;

            /**
             * Generate the local rhs vector
             * @param loc_vec - output vector
             * @param rhs_func - pointer to an rhs function
             * @param param - parameters of the task, some of them are used in rhs function
             * @param mesh_vertices - all mesh vertices
             * @param time - particular time for which we calculate this rhs vector
             */
  void local_rhs_vector(double *loc_vec,
                        double(*rhs_func)(const Point &point, double t, const Parameters &par),
                        const std::vector<Point> &points, double time, const Parameters &param) const;

            /**
             * Set the number of the DG edge
             * @param serial_num - serial number of the edge inside the triangle structure
             * @param num_edge - the number of edge itself
             */
  void con_edge(unsigned int serial_num, unsigned int num_edge);

            /**
             * Set the number of the DG edge
             * @param serial_num - serial number of the edge inside the triangle structure
             * @param num_edge - the number of edge itself
             */
  void dis_edge(unsigned int serial_num, unsigned int num_edge);

            /**
             * Check whether this triangle has this specific dof
             * @param num - the number of the degree of freedom
             */
  bool contains_dof(unsigned int num) const;


private: //========================== PRIVATE ============================

            /**
             * The number of partitions which the triangles is connected with.
             */
  std::vector<unsigned int> _ghost_cells;

            /**
             * Degrees of freedom defined on this triangle
             */
  std::vector<unsigned int> _dofs;

            /**
             * Determinant of the special matrix of coordinates of the triangles vertices
             */
  double _detD;

            /**
             * Coefficients _A of the barycentric functions on the triangle:
             * L[i] = _A[i]*x + _B[i]*y + _C[i]
             */
  double _A[n_vertices];

            /**
             * Coefficients _B of the barycentric functions on the triangle:
             * L[i] = _A[i]*x + _B[i]*y + _C[i]
             */
  double _B[n_vertices];

            /**
             * Coefficients _C of the barycentric functions on the triangle:
             * L[i] = _A[i]*x + _B[i]*y + _C[i]
             */
  double _C[n_vertices];

            /**
             * Edges for DG method. vector _edges from MeshElement class is used for CG edges
             */
  std::vector<unsigned int> _dis_edges;

            /**
             * Calculate all coefficients of barycentric functions and detD
             * @param mesh_vertices - the vertices of the mesh
             */
  void calculate_barycentric(const std::vector<Point> &mesh_vertices);

            /**
             * Get the gradient of one basis function specified by serial number (in triangle structure)
             * @param num_bf - serial number of the basis function of interest
             */
  const Point gradient_basis_function(unsigned int num_bf) const;
};




#endif // TRIANGLE_H
