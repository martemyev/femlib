#ifndef EDGE_H
#define EDGE_H

#include "config.h"
#include "line.h"
#include <vector>

NAMESPACE_FEM_OPEN


class Point;
class Triangle;


class Edge : public Line
{
public:
  Edge();
  ~Edge();

  Edge(const Edge &edge);
  const Edge& operator =(const Edge &edge);


            /**
             * Get the number of triangle which this edge belongs to
             */
  int triangle() const;

            /**
             * Set the number of triangle which this edge belongs to
             */
  void triangle(int num);

            /**
             * Get the constant reference to a list of edges that have the same physical
             * location that this edge does, but the different meaning.
             */
  const std::vector<unsigned int>& edges() const;

            /**
             * Set the list of edges that have the same physical
             * location that this edge does, but the different meaning.
             */
  void edges(const std::vector<unsigned int> &e);

  unsigned int n_assoc_edges() const;
  unsigned int assoc_edge(unsigned int num) const;

            /**
             * Calculate the length of the edges
             * @param points - a list of all mesh vertices (or, maybe, degrees of freedom)
             */
  double length(const std::vector<Point> &points) const;

            /**
             * Change the orientarion of the edge swaping its vertices
             */
  void swap_vertices();

            /**
             * Get an outward unit normal vector of the edge
             * @param tri - triangle which this edge belongs to. we need it to define outward direction of the vector
             * @param points - to know the coordinates of all vertices
             */
  const Point normal(const Triangle &tri, const std::vector<Point> &points) const;

  virtual double mes() const;


private:
            /**
             * The number of triangle which the edge belongs to.
             * This number exists for DG method only, since in this case
             * each triangle has its own edges. In the case of CG method
             * this field is not used, and filled by -1 by default.
             */
  int _triangle;

            /**
             * The numbers of edges that exists in DG method and coincide physically with
             * the edge from CG method. Therefore this field is valid for CG method only.
             * For DG method this field is not used.
             */
  std::vector<unsigned int> _edges;
};


NAMESPACE_FEM_CLOSE

#endif // EDGE_H
