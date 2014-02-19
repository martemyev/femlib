#include "edge.h"
#include "auxiliary_functions.h"
#include "point.h"
#include "triangle.h"
#include "math_functions.h"
#include <algorithm>


Edge::Edge()
  : _triangle(-1)
{ }



Edge::~Edge()
{
  _edges.clear();
}



Edge::Edge(const Edge &e)
  : Line(e),
    _triangle(e._triangle),
    _edges(e._edges)
{ }



const Edge& Edge::operator =(const Edge &e)
{
  Line::operator =(e);
  _triangle = e._triangle;
  _edges = e._edges;
  return *this;
}



void Edge::triangle(int num)
{
  expect(_triangle == -1, "The value of triangle is not default, that means it was already changed");
  _triangle = num;
}



void Edge::edges(const std::vector<unsigned int> &e)
{
  expect(e.size() == 1 || e.size() == 2, "Wrong size of the list of associated edges");
  _edges = e;
}



unsigned int Edge::n_assoc_edges() const
{
  return _edges.size();
}



unsigned int Edge::assoc_edge(unsigned int num) const
{
  expect(num < _edges.size(), "");
  return _edges[num];
}



int Edge::triangle() const
{
  return _triangle;
}



double Edge::length(const std::vector<Point> &points) const
{
  expect(_vertices[0] < points.size() && _vertices[1] < points.size(),
         "Given array of points doesn't suit to this edge, or the vertices of the edge are wrong");
  const double x0 = points[_vertices[0]].coord(0); // x-coordinate of the begin of the edge
  const double y0 = points[_vertices[0]].coord(1); // y-coordinate of the begin of the edge
  const double x1 = points[_vertices[1]].coord(0); // x-coordinate of the end of the edge
  const double y1 = points[_vertices[1]].coord(1); // y-coordinate of the end of the edge
  return sqrt((x0 - x1)*(x0 - x1) + (y0 - y1) * (y0 - y1));
}



void Edge::swap_vertices()
{
  expect(_vertices[0] != 0 || _vertices[1] != 0, "Vertices are not initialized");
  std::swap(_vertices[0], _vertices[1]);
}



const Point Edge::normal(const Triangle &tri, const std::vector<Point> &points) const
{
  expect(tri.contains_dof(_vertices[0]), "The first vertex (dof) of the edge is not in the triangle");
  expect(tri.contains_dof(_vertices[1]), "The second vertex (dof) of the edge is not in the triangle");

  const Point A = points[_vertices[0]]; // first vertex
  const Point B = points[_vertices[1]]; // second vertex
  // take the number of the third vertex
  int third_vertex = -1;
  for (int i = 0; i < Triangle::n_dofs_first && third_vertex == -1; ++i)
    if (tri.dof(i) != _vertices[0] &&
        tri.dof(i) != _vertices[1])
      third_vertex = tri.dof(i);
  expect(third_vertex != -1, "");
  const Point C = points[third_vertex]; // third vertex
  const Point AB = B - A; // vector connecting A and B
  const Point AC = C - A; // vector connecting A and C

  Point normal = Point(AB.coord(1), -AB.coord(0), 0); // normal vector has coordinates (AB.y, -AB.x);
  // check the direction of the normal
  if (dot_product(normal, AC) > 0) // if normal vector and AC have the same direction, then normal is oriented into triangle
    normal = Point(-AB.coord(1), AB.coord(0), 0); // in this case we just change the direction on outward

  // normalize the outward normal vector
  normalize(normal);

  return normal;
}
