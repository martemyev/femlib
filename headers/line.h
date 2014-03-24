#ifndef LINE_H
#define LINE_H

#include "config.h"
#include <vector>
#include "mesh_element.h"

NAMESPACE_FEM_OPEN



/**
 * Line keeps 3 numbers - the number of beginning vertex,
 * the number of the ending one, and a number of physical domain
 * where this line takes place (material identificator - in other words).
 * Line is not oriented.
 */
class Line : public MeshElement
{
public:
            /**
             * There are 2 vertices to describe a line
             */
  static const unsigned int n_vertices = 2;

            /**
             * Line is edge itself, so the number of edges is 1
             */
  static const unsigned int n_edges = 1;

            /**
             * It's 1D shape, so there is no faces here
             */
  static const unsigned int n_faces = 0;

            /**
             * In Gmsh line (physical line) is defined by number 1
             */
  static const unsigned int gmsh_el_type = 1;

            /**
             * Constructor
             */
  Line();

            /**
             * Constructor with parameters
             * @param ver - the list of vertices
             * @param mat_id - the material ID
             */
  Line(const std::vector<unsigned int> &ver,
       const unsigned int mat_id = 0);

            /**
             * Constructor with parameters
             * @param v1 - one vertex
             * @param v2 - another vertex
             * @param mat_id - material ID
             */
  Line(const unsigned int v1,
       const unsigned int v2,
       const unsigned int mat_id = 0);

            /**
             * Destructor
             */
  virtual ~Line();

            /**
             * Find common vertex between two lines
             * @param line - second line for seeking common vertex
             */
  unsigned int common_vertex(const Line& line) const;

            /**
             * Get another vertex (different from that we have)
             * @param vertex - we have the number of one vertex (this one),
             *                 and we want to find the number of another vertex
             */
  unsigned int another_vertex(const unsigned int vertex) const;

  virtual double mes() const;
};


NAMESPACE_FEM_CLOSE

#endif // LINE_H
