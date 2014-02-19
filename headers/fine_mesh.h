#ifndef FINE_MESH_H
#define FINE_MESH_H

#include "point.h"
#include "triangle.h"
#include "line.h"
#include <vector>
#include <string>
#include <set>

class Parameters;


class FineMesh
{
public:
            /**
             * Constructor
             */
  FineMesh(Parameters *param);

            /**
             * Destructor
             */
  ~FineMesh();

            /**
             * Read the mesh from a file
             * @param filename - a name of the file with the mesh
             */
  void read(const std::string &filename);

            /**
             * Get the number of mesh vertices
             */
  unsigned int n_vertices() const;
  void n_vertices(unsigned int amount);

            /**
             * Get the mesh vertex (its copy)
             * @param number - the serial number of the mesh vertex
             */
  Point vertex(unsigned int number) const;
  void vertex(unsigned int number, const Point &ver);

            /**
             * Get a constant reference to all mesh vertices
             */
  const std::vector<Point>& vertices() const;

            /**
             * Get the number of mesh triangles
             */
  unsigned int n_triangles() const;
  void n_triangles(unsigned int amount);

            /**
             * Get the mesh triangle (its copy)
             * @param number - the serial number of the mesh triangle
             */
  Triangle triangle(int number) const;
  void triangle(unsigned int number, const Triangle &tri);

            /**
             * Get the mesh triangle - not a copy - pointer to original one
             * @param number - the serial number of the mesh triangle
             */
  Triangle* triangle_orig(unsigned int number);

            /**
             * Get the maximal coordinates of the domain
             */
  Point max_coord() const;

            /**
             * Get the minimal coordinates of the domain
             */
  Point min_coord() const;

            /**
             * Get the number of the boundary lines
             */
  unsigned int n_lines() const;

            /**
             * Get the copy of the line
             * @param number - serial number of the line
             */
  Line line(unsigned int number) const;

            /**
             * Get the list of the boundary vertices of the mesh
             */
  const std::vector<int>& boundary_vertices() const;

            /**
             * Get the number of partitions of the mesh.
             * The partitions are associated (are considered) as coarse elements.
             */
  unsigned int n_partitions() const;

            /**
             * Check if the mesh has no vitally important elements (like nodes, triangles, etc)
             */
  bool empty() const;

            /**
             * Get the number of boundary vertices
             */
  unsigned int n_boundary_vertices() const;

            /**
             * Get a boundary vertex (i.e. the number of ordinary vertex which is boundary)
             * @param num - serial number of boundary vertex in the list
             */
  int boundary_vertex(unsigned int num) const;

  void boundary_vertices(const std::set<int> b_nodes);


private: // ========================== PRIVATE =========================
            /**
             * Parameters of the solving problem
             */
  Parameters *_param;

            /**
             * Mesh vertices
             */
  std::vector<Point> _vertices;

            /**
             * Maximal coordinates of the mesh vertices.
             * max_coord can not be one of the mesh vertices, if
             * the computational domain has curvilinear boundaries
             */
  Point _max_coord;

            /**
             * Minimal coordinates of the mesh vertices.
             * min_coord can not be one of the mesh vertices, if
             * the computational domain has curvilinear boundaries
             */
  Point _min_coord;

            /**
             * Mesh triangles
             */
  std::vector<Triangle> _triangles;

            /**
             * Mesh lines (like boundary lines)
             */
  std::vector<Line> _lines;

            /**
             * Mesh edges
             */
  std::vector<Line> _edges;

            /**
             * The numbers of the mesh vertices that are on the boundary of the computational domain
             */
  std::vector<int> _boundary_vertices;

            /**
             * The unique numbers of partitions of the mesh
             */
  std::set<unsigned int> _partitions;



  FineMesh(const FineMesh&);
  FineMesh& operator=(const FineMesh&);

            /**
             * Free the memory
             */
  void clear();

            /**
             * Generate the list of boundary vertices
             */
  void boundary_vertices_initialization();
};

#endif // FINE_MESH_H
