#ifndef MESH_ELEMENT_H
#define MESH_ELEMENT_H

#include <vector>


/**
 * This class implements the most part of functionality of
 * all elements of mesh: triangles, tetrahedra, quadrangles, hexahedra, etc.
 * All these elements are declared as pointers to base (this) class.
 * It's not an abstract class, because it has no pure virtual functions,
 * but you can't create objects of this class, because its constructor is protected.
 */
class MeshElement
{
public:

            /**
             * Destructor
             */
  virtual ~MeshElement();

            /**
             * Get the number of vertices
             */
  unsigned int n_vertices() const;

            /**
             * Get the number of edges
             */
  unsigned int n_edges() const;

            /**
             * Get the number of faces
             */
  unsigned int n_faces() const;

            /**
             * Get type of the element that is used in Gmsh
             */
  unsigned int gmsh_el_type() const;

            /**
             * Get the material ID of the element.
             * It's a number that describes the physical domain
             * which the element belongs to.
             */
  unsigned int material_id() const;

            /**
             * Get the number of partition which the element belongs to.
             * Partitions can be associated with coarse elements.
             */
  unsigned int partition_id() const;

            /**
             * Get the number of vertex describing the element
             * @param number - local number of vertex [0, n_vertices)
             * @return global number of vertex (among other mesh vertices)
             */
  unsigned int vertex(unsigned int number) const;

            /**
             * Get the number of edge describing the element
             * @param number - local number of edge [0, n_edges)
             * @return global number of edge (among other mesh edges)
             */
  unsigned int edge(unsigned int number) const;

            /**
             * Get the number of face describing the element
             * @param number - local number of face [0, n_faces)
             * @return global number of face (among other mesh faces)
             */
  unsigned int face(unsigned int number) const;

            /**
             * Set the number of vertex
             * @param local_number - the number of vertex inside the element [0, n_vertices)
             * @param global_number - the number of vertex among other vertices of the mesh
             */
  void vertex(unsigned int local_number, unsigned int global_number);

            /**
             * Set the number of edge
             * @param local_number - the number of edge inside the element [0, n_edges)
             * @param global_number - the number of edge among other edges of the mesh
             */
  void edge(unsigned int local_number, unsigned int global_number);

            /**
             * Set the number of face
             * @param local_number - the number of face inside the element [0, n_faces)
             * @param global_number - the number of face among other faces of the mesh
             */
  void face(unsigned int local_number, unsigned int global_number);

            /**
             * Set all faces once at time
             * @param face_numbers - the numbers of all cell faces
             */
  void faces(const std::vector<unsigned int> &face_numbers);

            /**
             * Check - whether the element contains the vertex or not
             * @param vertex - the number of vertex that we want to check
             */
  bool contains(const unsigned int vertex) const;


protected:
            /**
             * The number of vertices describing the element.
             * It must be defined in each derived class,
             * because it's 0 by default.
             */
  unsigned int _n_vertices;

            /**
             * Vertices (i.e. their global numbers) describing the element
             */
  std::vector<unsigned int> _vertices;

            /**
             * The number of edges describing the element.
             * It must be defined in each derived class,
             * because it's 0 by default.
             */
  unsigned int _n_edges;

            /**
             * Edges (i.e. their global numbers) describing the element
             * It's not always used.
             */
  std::vector<unsigned int> _edges;

            /**
             * The number of faces describing the element.
             * It must be defined in each derived class,
             * because it's 0 by default.
             */
  unsigned int _n_faces;

            /**
             * Faces (i.e. their global numbers) describing the element
             * It's not always used.
             */
  std::vector<unsigned int> _faces;

            /**
             * ID of the physical domain where the element takes place.
             * It's necessary to distinguish media with different physical properties.
             */
  unsigned int _material_id;

            /**
             * ID of partition (like domain decomposition)
             * which can be considered as a coarse element
             */
  unsigned int _partition_id;

            /**
             * Type of the element (its number actually) like in Gmsh.
             * It must be defined in every derived class.
             * It's 0 by default.
             */
  unsigned int _gmsh_el_type;

            /**
             * Constructor is protected to prevent creating MeshElement objects directly
             * @param n_ver - number of vertices
             * @param n_edg - number of edges
             * @param n_fac - number of faces
             * @param el_type - type of the element in Gmsh
             */
  MeshElement(unsigned int n_ver = 0,
              unsigned int n_edg = 0,
              unsigned int n_fac = 0,
              unsigned int el_type = 0);

            /**
             * Copy constructor
             */
  MeshElement(const MeshElement &elem);

            /**
             * Copy assignment operator
             */
  MeshElement& operator =(const MeshElement &elem);
};


#endif // MESH_ELEMENT_H
