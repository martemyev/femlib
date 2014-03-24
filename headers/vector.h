#ifndef FEM_VECTOR_H
#define FEM_VECTOR_H

#include "config.h"
#include "petscvec.h"
#include <vector>

NAMESPACE_FEM_OPEN



/**
 * Wrapper for PETSc Vec class
 */
class Vector
{
public:
  Vector();
  Vector(const Vector &vec);
  Vector& operator =(const Vector &vec);

  Vector(int size);

  ~Vector();

            /**
             * Set the size of the vector
             * @param s - the size
             */
  void size(int s);

            /**
             * Get the size of the vector
             */
  int size() const;

            /**
             * Get the i-th element
             * @param i - the number of the element we want to return
             */
  double operator ()(int i) const;

            /**
             * Extract all values of the vector. It can be also considered as converting Vec class to the std::vector
             */
  void values(std::vector<double> &v) const;

  void add(int i, double value);
  void insert(int i, double value);

  const Vec& vec() const;
  void init(const Vec &vector);


private:
  int _size;
  Vec _data;
};


NAMESPACE_FEM_CLOSE

#endif // FEM_VECTOR_H
