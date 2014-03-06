#ifndef FEM_RECTANGLE_H
#define FEM_RECTANGLE_H

#include "quadrangle.h"


class Function;


/**
 * Rectangle is a quadrangle that has right angles (Pi/2)
 */
class Rectangle : public Quadrangle
{
public:
            /**
             * Default constructor
             */
  Rectangle();

            /**
             * Virtual destructor.
             * It does nothing but passes the action to the Quadrangle's destructor.
             */
  virtual ~Rectangle();

            /**
             * Constructor with parameters
             * @param ver - rectangle's vertices
             * @param mesh_vertices - all vertices of the mesh.
             *                        they are needed to save the coordinates of the vertices
             *                        for faster computation of basis functions in the future
             * @param mat_id - material ID
             * @param part_id - partition ID (can be used as a number of a coarse element)
             * @param ghost_cells - the numbers of other partitions (coarse elements)
             *                      that are connected with this rectangle at least by one vertex
             */
  Rectangle(const std::vector<unsigned int> &ver,
            const std::vector<Point> &mesh_vertices = std::vector<Point>(),
            const unsigned int mat_id = 0,
            const unsigned int part_id = 0,
            const std::vector<unsigned int> &ghost_cells = std::vector<unsigned int>());

            /**
             * Constructor with parameters
             * @param v1 - first vertex
             * @param v2 - second vertex
             * @param v3 - third vertex
             * @param v4 - fourth vertex
             * @param mesh_vertices - all vertices of the mesh.
             *                        they are needed to save the coordinates of the vertices
             *                        for faster computation of basis functions in the future
             * @param mat_id - material ID
             * @param part_id - partition ID (can be used as a number of a coarse element)
             * @param ghost_cells - the numbers of other partitions (coarse elements)
             *                      that are connected with this rectangle at least by one vertex
             */
  Rectangle(const unsigned int v1,
            const unsigned int v2,
            const unsigned int v3,
            const unsigned int v4,
            const std::vector<Point> &mesh_vertices = std::vector<Point>(),
            const unsigned int mat_id = 0,
            const unsigned int part_id = 0,
            const std::vector<unsigned int> &ghost_cells = std::vector<unsigned int>());

            /**
             * Copy constructor
             */
  Rectangle(const Rectangle &rect);

            /**
             * Copy assignment operator
             */
  Rectangle& operator =(const Rectangle &rect);

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
             * Generate the local rhs vector
             * @param rhs_func - pointer to an rhs function
             * @param points - all mesh vertices (or relevant points, like degrees of freedom)
             * @param time - particular time for which we calculate this rhs vector
             * @param loc_vec - output vector
             */
  void local_rhs_vector(const Function &func,
                        const std::vector<Point> &points,
                        const double time,
                        double *loc_vec) const;

  unsigned int n_dofs() const;
  void n_dofs(unsigned int number);
  unsigned int dof(unsigned int number) const;
  void dof(unsigned int number, unsigned int value);



private: // ========================= PRIVATE ====================
            /**
             * Degrees of freedom defined on this rectangle
             */
  std::vector<unsigned int> _dofs;

            /**
             * Check that this is rectangle and not just quadrangle,
             * i.e. that the angles are right
             */
  void check_rectangle() const;
};


#endif // FEM_RECTANGLE_H
