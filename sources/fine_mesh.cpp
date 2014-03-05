#include "fine_mesh.h"
#include "auxiliary_functions.h"
#include "math_functions.h"
#include <fstream>
#include <map>
#include <algorithm>


FineMesh::FineMesh()
{ }



FineMesh::~FineMesh()
{
  clear();
}



void FineMesh::clear()
{
  _vertices.clear();
  _triangles.clear();
  _lines.clear();
  _edges.clear();
  _boundary_vertices.clear();
  _partitions.clear();
}



void FineMesh::read(const std::string &file,
                    const Point &declared_min_point,
                    const Point &declared_max_point)
{
  std::ifstream in(file.c_str());
  require(in, "File " + file + " cannot be opened!");

  clear(); // if some mesh has been already read before, we delete it

  std::string str;
  in >> str; // the first string of Gmsh file is "$MeshFormat"
  expect(str == "$MeshFormat",
         "The first string of the Gmsh file " + file +
         " doesn't equal to \"$MeshFormat\". The actual string is \"" + str + "\"");

  // read the information about the mesh
  double version;
  int binary, dsize;
  in >> version >> binary >> dsize;
  expect(version == 2.2,
         "The version of Gmsh's mesh is unexpected (" + d2s(version) + ").");
  expect(dsize == sizeof(double),
         "The size of Gmsh's double (" + d2s(dsize) +\
         ") doesn't equal to size of double type (" + d2s(sizeof(double)) + ")");

  getline(in, str); // read some empty string

  // there is additional 1 (the number - one) in binary format
  if (binary)
  {
    int one;
    in.read(reinterpret_cast<char*>(&one), sizeof(int));
    require(one == 1,
            "The binary one (" + d2s(one) + ") doesn't equal to 1!");
  }

  // we make a map between serial number of the vertex and its number in the file.
  // it will help us when we create mesh elements
  std::map<unsigned int, unsigned int> vertices_map;

  // read lines of mesh file.
  // if we face specific keyword, we'll treat the section.
  while (in >> str)
  {
    if (str == "$Nodes") // read the mesh vertices
    {
      unsigned int n_vertices; // the number of all mesh vertices (that are saved in the file)
      in >> n_vertices; // read that number
      _vertices.resize(n_vertices); // allocate the memory for mesh vertices
      getline(in, str); // read some empty string

      unsigned int number; // the number of the vertex
      double coord[Point::n_coord]; // Cartesian coordinates of the vertex (Gmsh produces 3D mesh regardless its real dimension)
      double max_coord[Point::n_coord]; // the limits of the computational domain: maximum
      double min_coord[Point::n_coord]; // and minimum

      // read vertices
      for (unsigned int ver = 0; ver < n_vertices; ++ver)
      {
        if (binary) // binary format
        {
          in.read(reinterpret_cast<char*>(&number), sizeof(unsigned int)); // the number of each node
          in.read(reinterpret_cast<char*>(coord), Point::n_coord * sizeof(double)); // node coordinates
        }
        else // ASCII format
        {
          in >> number;
          for (unsigned int i = 0; i < Point::n_coord; ++i)
            in >> coord[i];
          if (ver == 0) // for the first vertex
          {
            for (unsigned int i = 0; i < Point::n_coord; ++i)
              max_coord[i] = min_coord[i] = coord[i]; // initialization of the max and min coordinates
          }
          else // for the other vertices
          {
            for (unsigned int i = 0; i < Point::n_coord; ++i)
            {
              max_coord[i] = std::max(max_coord[i], coord[i]); // searching max and min coordinates
              min_coord[i] = std::min(min_coord[i], coord[i]);
            }
          }

        }
        _vertices[ver] = Point(coord); // save the vertex
        vertices_map[number] = ver; // add the number of vertex to the map
      }
      _max_coord = Point(max_coord); // this point can not be one of the mesh vertices if the domain has curvilinear boundaries
      _min_coord = Point(min_coord); // the same as above is true

      expect(n_vertices == vertices_map.size(),
             "Vertices numbers are not unique: n_vertices = " + d2s(n_vertices) +
             " vertices_map.size() = " + d2s(vertices_map.size()));

    } // read the vertices

    else if (str == "$Elements") // read the mesh elements
    {
      unsigned int n_elements; // the number of mesh elements
      in >> n_elements; // read that number
      getline(in, str); // empty string

      unsigned int number; // the number of the element [1, nElements]
      unsigned int el_type; // the type of the element (1 - line, 2 - triangle, etc)
      unsigned int n_tags; // the number of tags describing the element
      unsigned int phys_domain; // the physical domain where the element takes place
      unsigned int elem_domain; // the elementary domain where the element takes place
      unsigned int n_partitions; // the number of partitions in which the element takes place
      unsigned int partition; // the partition which the element belongs to
      std::vector<unsigned int> ghost_cells; // "ghost cells" are other partitions which this element is connected with

      // the map between the type of the element,
      // and the number of nodes that describe it
      std::map<unsigned int, unsigned int> type_nodes;
      type_nodes[1] = 2; // 2-nodes line
      type_nodes[2] = 3; // 3-nodes triangle
      type_nodes[3] = 4; // 4-nodes quadrangle
      type_nodes[4] = 4; // 4-nodes tetrahedron
      type_nodes[5] = 8; // 8-nodes hexahedron
      type_nodes[15]= 1; // 1-node point

      if (binary) // binary format
      {
        require(false, "Binary 2.0 format of Gmsh mesh files is not supported");
      } // binary format

      else // ASCII format
      {
        for (int el = 0; el < n_elements; ++el)
        {
          in >> number >> el_type >> n_tags;
          std::vector<int> data(n_tags); // allocate the memory for some data
          for (unsigned int i = 0; i < n_tags; ++i) // read this information
            in >> data[i];
          phys_domain = (n_tags > 0) ? data[0] : 0; // physical domain - the most important value
          elem_domain = (n_tags > 1) ? data[1] : 0; // elementary domain
          if (n_tags > 2)
          {
            // the number of partitions where this elements takes place
            n_partitions = data[2];
            expect(n_partitions >= 1,
                   "The number of tags is more than 2. That means that we have partitions. But the number of partitions is "
                   + d2s(n_partitions));
            // the partition which the element belongs to
            partition = data[3] - 1; // since we associate the number of partition (which is numerated from 1) with the number of coarse element (which is numerated from 0)
            _partitions.insert(partition);
            // "ghost cells"
            if (n_partitions > 1) // if the element is on the boundary between the partitions, it is described by "ghost cells" as well
            {
              ghost_cells.resize(n_partitions - 1);
              for (int gc = 0; gc < n_partitions - 1; ++gc)
              {
                ghost_cells[gc] = -data[4 + gc]; // 'minus' since ghost cells are described by number of partition with the negative sing
                expect(ghost_cells[gc] > 0, "The number of the ghost cell (positive one) is unexpected (" + d2s(ghost_cells[gc]) + ")");
                --ghost_cells[gc]; // the same reason as in case of the number of partition
              }
            }
          }

          data.clear(); // clear the memory

          // how many vertices (nodes) describe the element
          std::map<unsigned int, unsigned int>::const_iterator el_type_iter =
              type_nodes.find(el_type);

          require(el_type_iter != type_nodes.end(),
                  "This type of the Gmsh's element (" + d2s(el_type) +
                  ") in the mesh file \"" + file + "\" is unknown!");

          const unsigned int n_elem_nodes = el_type_iter->second; // the number of nodes
          std::vector<unsigned int> nodes(n_elem_nodes); // allocate memory for nodes
          for (unsigned int i = 0; i < n_elem_nodes; ++i)
          {
            in >> nodes[i]; // read the numbers of nodes
            // vertices can be numerated not sequentially (or at least not from 0)
            nodes[i] = vertices_map.find(nodes[i])->second;
          }

          // add new element in the list
          //MeshElement *new_element;
          switch (el_type)
          {
//          case 15: // 1-node point
//            points.push_back(MeshElement_ptr(new PhysPoint(nodes, phys_domain)));
//            break;
          case 1: // 2-nodes line
            _lines.push_back(Line(nodes, phys_domain));
            //new_element = new Line(nodes, phys_domain);
            break;
          case 2: // 3-nodes triangle
            _triangles.push_back(Triangle(nodes, _vertices, phys_domain, partition, ghost_cells));
            //new_element = new Triangle(nodes, phys_domain);
            break;
//          case 3: // 4-nodes quadrangle
//            quadrangles.push_back(MeshElement_ptr(new Quadrangle(nodes, phys_domain)));
//            //new_element = new Quadrangle(nodes, phys_domain);
//            break;
//          case 4: //4-nodes tetrahedron
//            tetrahedra.push_back(MeshElement_ptr(new Tetrahedron(nodes, phys_domain)));
//            //new_element = new Tetrahedron(nodes, phys_domain);
//            break;
//          case 5: // 8-nodes hexahedron
//            hexahedra.push_back(MeshElement_ptr(new Hexahedron(nodes, phys_domain)));
//            //new_element = new Hexahedron(nodes, phys_domain);
//            break;
          default:
            require(false,
                    "Unknown type of the Gmsh's element (" + d2s(el_type) +
                    ") in the file " + file + "!");
          }

          nodes.clear();

          //elements.push_back(new_element);
          //delete new_element;
        }

        // check some expectations
        expect(number == n_elements,
               "The number of the last read Gmsh's element (" + d2s(number) +\
               ") is not equal to the amount of all elements in the mesh (" + d2s(n_elements) + ")!");

      } // ASCII format

      // requirements after reading elements
      require(!_triangles.empty(),
             "There are no any 2D or 3D elements in the mesh!");
      //require(!_lines.empty(),
      //       "There are no boundary lines in the mesh!");

      // check that the numbers of partitions are sequential
      for (int par = 0; par < _partitions.size(); ++par)
        require(*_partitions.find(par) == par, "The numeration of the partitions is not dense");


    } // read the elements
  } // read the mesh file

  in.close(); // close the file

  // generate the list of boundary vertices
  boundary_vertices_initialization(declared_min_point, declared_max_point);
}



void FineMesh::boundary_vertices_initialization(const Point &declared_min_point,
                                                const Point &declared_max_point)
{
  if (_lines.size() == 0) // there are no boundary lines in the mesh
  {
    require(norm(declared_max_point - declared_min_point) > FLOAT_NUMBERS_EQUALITY_TOLERANCE,
            "There are no boundary lines in the mesh, and boundaries seem not to be defined through min and max points");
    // in this case the domain should be rectangular.
    const double tol = FLOAT_NUMBERS_EQUALITY_TOLERANCE;
    require(fabs(_min_coord.coord(0) - declared_min_point.coord(0)) < tol &&
            fabs(_min_coord.coord(1) - declared_min_point.coord(1)) < tol &&
            fabs(_max_coord.coord(0) - declared_max_point.coord(0)) < tol &&
            fabs(_max_coord.coord(1) - declared_max_point.coord(1)) < tol,
            "There are no boundary lines and declared computational domain: (x0,y0)x(x1,y1)=(" +
            d2s(declared_min_point.coord(0)) + ", " + d2s(declared_min_point.coord(1)) + ")x(" +
            d2s(declared_max_point.coord(0)) + ", " + d2s(declared_max_point.coord(1)) +
            ") is not equal to the domain read from the mesh file: (" +
            d2s(_min_coord.coord(0)) + ", " + d2s(_min_coord.coord(1)) + ")x(" +
            d2s(_max_coord.coord(0)) + ", " + d2s(_max_coord.coord(1)) + ")");

    // find nodes that lie on the boundary of the computational domain
    for (int i = 0; i < _vertices.size(); ++i)
    {
      const double x = _vertices[i].coord(0);
      const double y = _vertices[i].coord(1);
      if ((fabs(x - _min_coord.coord(0)) < tol ||
           fabs(y - _min_coord.coord(1)) < tol ||
           fabs(x - _max_coord.coord(0)) < tol ||
           fabs(y - _max_coord.coord(1)) < tol) &&
          find(_boundary_vertices.begin(), _boundary_vertices.end(), i) == _boundary_vertices.end())
        _boundary_vertices.push_back(i);
    }
  }
  else // there are the boundary lines in the mesh
  {
    for (int lin = 0; lin < _lines.size(); ++lin)
    {
      for (int ver = 0; ver < Line::n_vertices; ++ver)
      {
        const int node = _lines[lin].vertex(ver);
        if (find(_boundary_vertices.begin(), _boundary_vertices.end(), node) == _boundary_vertices.end())
          _boundary_vertices.push_back(node);
      }
    }
  }
}



unsigned int FineMesh::n_vertices() const
{
  return _vertices.size();
}



void FineMesh::n_vertices(unsigned int amount)
{
  expect(_vertices.empty(), "The list of vertices is not empty. What are you doing?");
  _vertices.resize(amount); // allocate the memory for the vector
}



Point FineMesh::vertex(unsigned int number) const
{
  expect(number >= 0 && number < _vertices.size(), "Incorrect input parameter");
  return _vertices[number];
}



void FineMesh::vertex(unsigned int number, const Point &ver)
{
  expect(number < _vertices.size(), "Incorrect input parameter");
  _vertices[number] = ver;
}



const std::vector<Point>& FineMesh::vertices() const
{
  return _vertices;
}



unsigned int FineMesh::n_triangles() const
{
  return _triangles.size();
}



void FineMesh::n_triangles(unsigned int amount)
{
  expect(_triangles.empty(), "The list of triangles is not empty. What are you doing?");
  _triangles.resize(amount); // allocate the memory for the vector
}



Triangle FineMesh::triangle(int number) const
{
  expect(number >= 0 && number < _triangles.size(),
         "Incorrect number of a triangle (" + d2s(number) +
         "). It should be in a range [0, " + d2s(_triangles.size()) + ")");
  return _triangles[number];
}



void FineMesh::triangle(unsigned int number, const Triangle &tri)
{
  expect(number < _triangles.size(),
         "Incorrect number of a triangle (" + d2s(number) +
         "). It should be in a range [0, " + d2s(_triangles.size()) + ")");
  _triangles[number] = tri;
}



Triangle* FineMesh::triangle_orig(unsigned int number)
{
  expect(number >= 0 && number < _triangles.size(), "Incorrect input parameter");
  return &_triangles[number];
}



Point FineMesh::max_coord() const
{
  return _max_coord;
}



Point FineMesh::min_coord() const
{
  return _min_coord;
}



unsigned int FineMesh::n_lines() const
{
  return _lines.size();
}



Line FineMesh::line(unsigned int number) const
{
  expect(number >= 0 && number < _lines.size(), "Incorrect input parameter");
  return _lines[number];
}



std::vector<int> const& FineMesh::boundary_vertices() const
{
  expect(!_boundary_vertices.empty(), "The vector of boundary vertices is empty");
  return _boundary_vertices;
}



unsigned int FineMesh::n_partitions() const
{
  return _partitions.size();
}



bool FineMesh::empty() const
{
  return (_vertices.empty() && (_triangles.empty() || _rectangles.empty()));
}



unsigned int FineMesh::n_boundary_vertices() const
{
  return _boundary_vertices.size();
}



int FineMesh::boundary_vertex(unsigned int num) const
{
  expect(num < _boundary_vertices.size(), "Incorrect input parameter");
  return _boundary_vertices[num];
}



void FineMesh::boundary_vertices(const std::set<int> b_nodes)
{
  expect(_boundary_vertices.empty(), "The list of boundary vertices is not empty. WTF?");
  _boundary_vertices.resize(b_nodes.size()); // allocate the memory
  std::copy(b_nodes.begin(), b_nodes.end(), _boundary_vertices.begin()); // copy the values
}



const std::vector<Rectangle>& FineMesh::rectangles() const
{
  return _rectangles;
}



unsigned int FineMesh::n_rectangles() const
{
  return _rectangles.size();
}



Rectangle* FineMesh::rectangle_orig(unsigned int number)
{
  expect(number < _rectangles.size(), "Incorrect input parameter");
  return &_rectangles[number];
}



void FineMesh::create_rectangular_grid(double X_BEG, double X_END,
                                       double Y_BEG, double Y_END,
                                       unsigned int N_FINE_X, unsigned int N_FINE_Y)
{
  require(X_END > X_BEG && Y_END > Y_BEG, "");

  const double hx = (X_END - X_BEG) / (double)N_FINE_X; // step in x-direction
  const double hy = (Y_END - Y_BEG) / (double)N_FINE_Y; // step in y-direction

  const unsigned int n_vertices = (N_FINE_X + 1) * (N_FINE_Y + 1); // the total number of vertices
  _vertices.resize(n_vertices); // allocate the memory for all vertices

  const unsigned int n_rectangles = N_FINE_X * N_FINE_Y; // the total number of rectangles
  _rectangles.resize(n_rectangles); // allocate the memory for all rectangles

  unsigned int ver = 0; // number of a current vertex
  for (int i = 0; i < N_FINE_Y + 1; ++i)
  {
    const double y = (i == N_FINE_Y ? Y_END : Y_BEG + i * hy);
    for (int j = 0; j < N_FINE_X + 1; ++j)
    {
      const double x = (j == N_FINE_X ? X_END : X_BEG + j * hx);
      _vertices[ver] = Point(x, y, 0);
      ++ver;
    }
  }

  unsigned int rect = 0; // number of a current rectangle
  std::vector<unsigned int> vert(Rectangle::n_vertices); // the numbers of vertices describing a rectangle
  for (int i = 0; i < N_FINE_Y; ++i)
  {
    for (int j = 0; j < N_FINE_X; ++j)
    {
      // numeration of rectangle's vertices is the following
      // 2 --- 3
      // |     |
      // |     |
      // 0 --- 1
      vert[0] = i * (N_FINE_X + 1) + j;
      vert[1] = i * (N_FINE_X + 1) + j + 1;
      vert[2] = (i + 1) * (N_FINE_X + 1) + j;
      vert[3] = (i + 1) * (N_FINE_X + 1) + j + 1;
      _rectangles[rect] = Rectangle(vert, _vertices);
      ++rect;
    }
  }

  // boundary vertices initialization
  _min_coord = Point(X_BEG, Y_BEG);
  _max_coord = Point(X_END, Y_END);
  boundary_vertices_initialization(_min_coord, _max_coord);
}



Rectangle FineMesh::rectangle(unsigned int number) const
{
  expect(number < _rectangles.size(), "");
  return _rectangles[number];
}
