#ifndef FEM_QUADRANGLE_H
#define FEM_QUADRANGLE_H

#include "mesh_element.h"
#include "point.h"


/**
 * Quadrangle is a 2-dimensional shape with 4 straight edges.
 * It's an element of mesh,
 * therefore it inherits the most part of
 * functionality from MeshElement class.
 */
class Quadrangle : public MeshElement
{
public:
            /**
             * The number of vertices of quadrangle
             */
  static const unsigned int n_vertices = 4;

            /**
             * The number of edges of quadrangle
             */
  static const unsigned int n_edges = 4;

            /**
             * Quadrangle is a 2D shape,
             * so it's a face itself (for hexahedron)
             */
  static const unsigned int n_faces = 1;

            /**
             * In Gmsh a quadrangle is defined by number 3
             */
  static const unsigned int gmsh_el_type = 3;

            /**
             * In VTK format a triangle is defined by number 9 (VTK_QUAD)
             */
  static const unsigned int vtk_el_type = 9;

            /**
             * The number of degrees of freedom on the triangle
             * in case of first order basis functions coincides
             * with the number of vertices
             */
  static const unsigned int n_dofs_first = n_vertices;

            /**
             * Default constructor
             */
  Quadrangle();

            /**
             * Virtual destructor.
             * It does nothing but passes the action to the MeshElement's destructor.
             */
  virtual ~Quadrangle();

            /**
             * Constructor with parameters
             * @param ver - quadrangle vertices
             * @param mesh_vertices - all vertices of the mesh.
             *                        they are needed to save the coordinates of the vertices
             *                        for faster computation of basis functions in the future
             * @param mat_id - material ID
             * @param part_id - partition ID (can be used as a number of a coarse element)
             * @param ghost_cells - the numbers of other partitions (coarse elements)
             *                      that are connected with this quadrangle at least by one vertex
             */
  Quadrangle(const std::vector<unsigned int> &ver,
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
             *                      that are connected with this quadrangle at least by one vertex
             */
  Quadrangle(const unsigned int v1,
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
  Quadrangle(const Quadrangle &quad);

            /**
             * Copy assignment operator
             */
  Quadrangle& operator =(const Quadrangle &quad);


protected: // ======================= PROTECTED =======================
            /**
             * x-coordinates of the vertices of the quadrangle.
             * they are saved for faster computation in the future
             */
  double _X[n_vertices];

            /**
             * y-coordinates of the vertices of the quadrangle.
             * they are saved for faster computation in the future
             */
  double _Y[n_vertices];

};



#endif // FEM_QUADRANGLE_H
