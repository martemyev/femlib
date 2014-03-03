#ifndef FEM_QUADRANGLE_H
#define FEM_QUADRANGLE_H

#include "mesh_element.h"

/**
 * Quadrangle - 2-dimensional shape with 4 straight edges.
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
             * Quadrangle is 2D shape,
             * so it's a face itself (for hexahedron)
             */
  static const unsigned int n_faces = 1;

            /**
             * In Gmsh quadrangle is defined by number 3
             */
  static const unsigned int gmsh_el_type = 3;

            /**
             * Default constructor.
             */
  Quadrangle();

            /**
             * Constructor with parameters
             * @param ver - quadrangle vertices
             * @param mat_id - material ID
             */
  Quadrangle(const std::vector<unsigned int> &ver,
             const unsigned int mat_id = 0);

            /**
             * Constructor with parameters
             * @param v1 - first vertex
             * @param v2 - second vertex
             * @param v3 - third vertex
             * @param v4 - fourth vertex
             * @param mat_id - material ID
             */
  Quadrangle(const unsigned int v1,
             const unsigned int v2,
             const unsigned int v3,
             const unsigned int v4,
             const unsigned int mat_id = 0);

            /**
             * Virtual destructor.
             * It does nothing passing the action to MeshElement's desctructor.
             */
  virtual ~Quadrangle();

private: // ======================= PRIVATE =======================


};



#endif // FEM_QUADRANGLE_H
