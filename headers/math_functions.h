#ifndef MATH_FUNCTIONS_H
#define MATH_FUNCTIONS_H

#include "config.h"
#include "petscvec.h"
#include "point.h"

NAMESPACE_FEM_OPEN

namespace math {



const double FLOAT_NUMBERS_EQUALITY_TOLERANCE = 1e-14;
const double FLOAT_NUMBERS_EQUALITY_REDUCED_TOLERANCE = 1e-6;
const double PI = 3.141592654;


/**
 * Relative error in 2-norm between two vectors
 * rel_error = norm2(vec1-vec2) / norm2(vec2)
 */
double rel_error(const Vec &vec1, const Vec &vec2);


/**
 * Dot (inner) product of two 3D vectors
 * @param p1 - one vector (represented as object of class Point)
 * @param p2 - another vector (represented as object of class Point)
 */
double dot_product(const Point &p1, const Point &p2);


/**
 * Vector normalization. The vector is changed
 * @param vec - vector (represented as object of class Point)
 * @return
 */
void normalize(Point &vec);



/**
 * Norm of a vector represented by a Point object
 */
double norm(const Point &vec);


} // namespace math

NAMESPACE_FEM_CLOSE

#endif // MATH_FUNCTIONS_H
