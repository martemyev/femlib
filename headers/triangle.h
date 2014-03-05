#ifndef FEM_TRIANGLE_H
#define FEM_TRIANGLE_H

#include "mesh_element.h"
#include "point.h"


class Edge;
class DoFHandler;


/**
 * Triangle is a 2-dimensional simplex.
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
             * In VTK format a triangle is defined by number 5 (VTK_TRIANGLE)
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
             * @param mesh_vertices - all vertices of the mesh. they are needed for calculation of barycentric functions
             * @param mat_id - material ID
             * @param part_id - partition ID (can be used as a number of a coarse element)
             * @param ghost_cells - the numbers of other partitions (coarse elements)
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
             * @param mesh_vertices - all vertices of the mesh. they are needed for calculation of barycentric functions
             * @param mat_id - material ID
             * @param part_id - partition ID (can be used as a number of a coarse element)
             * @param ghost_cells - the numbers of other partitions (coarse elements)
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
             * Get the number of the ghost cells, i.e. the number of connected triangle from another partition (coarse element).
             * Connection means that the triangles have at least one common vertex.
             */
  unsigned int n_ghost_cells() const;

            /**
             * Get the number of 'num'-th ghost cell
             * @param num - the serial number of a ghost cell from the vector of them
             */
  unsigned int ghost_cell(unsigned int num) const;

            /**
             * Get the number of degrees of freedom (dofs) associated with this triangle.
             * That's the real number of dofs allocated for this particular triangle.
             * It can be that a triangle has no dofs at all. That's possible if the procedure of dofs distribution hasn't been passed yet.
             * There are some constants though (like n_dofs_first),
             * that define the number of degrees of freedom that a triangle
             * should have according to the order of basis functions defined on it.
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
             * Generate the local mass matrix for the triangle K
             * M_{ij}^K = \int_K \alpha \phi_i \phi_j dK,
             * where \phi are basis functions
             * @param coef_alpha - constant coefficient under an integral
             * @param loc_mat - output data, generated local matrix
             */
  void local_mass_matrix(const double coef_alpha, double **loc_mat) const;

            /**
             * Generate the local stiffness matrix for the triangle K
             * S_{ij}^K = \int_K \beta \nabla\phi_i \cdot \nabla\phi_j dK,
             * where \phi are basis functions
             * @param coef_beta - constant coefficient under an integral
             * @param loc_mat - output data, generated local matrix
             */
  void local_stiffness_matrix(const double coef_beta, double **loc_mat) const;

            /**
             * Generate the local matrix of integrals over a one boundary edge belogning to this triangle
             * @param edge - the edge over which we need to take the integrals
             * @param coef_beta - constant coefficient under an integral
             * @param gamma - the IPDG parameter
             * @param dof_handler - the constant reference to an object of DoFHandler class
             *                      that provides us with the list of all degrees of freedom, for example
             * @param loc_mat - output data, generated local matrix
             */
  void local_dg_boundary_matrix(const Edge &edge,
                                const double coef_beta,
                                const double gamma,
                                const DoFHandler &dof_handler,
                                double **loc_mat) const;

            /**
             * Generate the local matrix of integrals over a one interior CG edge sharing by this triangle and another one
             * @param tri - another triangle
             * @param edges - the DG edges: edges[0] belongs to this triangle, edges[1] - to another one
             * @param coefs - the coefficients beta (under an integral): coefs[0] is on this triangle, coefs[1] is on another one
             * @param gamma - the IPDG parameter
             * @param dof_handler - the constant reference to an object of DoFHandler class
             *                      that provides us with the list of all degrees of freedom, for example
             * @param loc_mat - output data, generated local matrix
             */
  void local_dg_interior_matrix(const Triangle &tri,
                                Edge edges[],
                                double coefs[],
                                const double gamma,
                                const DoFHandler &dof_handler,
                                double **loc_mat) const;

            /**
             * Generate the local rhs vector
             * @param rhs_func - pointer to an rhs function
             * @param points - all mesh vertices (or relevant points, like degrees of freedom)
             * @param time - particular time for which we calculate this rhs vector
             * @param loc_vec - output vector
             */
  void local_rhs_vector(double(*rhs_func)(const Point &point, double t),
                        const std::vector<Point> &points,
                        const double time,
                        double *loc_vec) const;

            /**
             * Set the number of the CG edge
             * @param serial_num - serial number of the edge inside the triangle's structure
             * @param num_edge - the number of the edge itself
             */
  void con_edge(unsigned int serial_num, unsigned int num_edge);

            /**
             * Set the number of the DG edge
             * @param serial_num - serial number of the edge inside the triangle's structure
             * @param num_edge - the number of the edge itself
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




#endif // FEM_TRIANGLE_H
